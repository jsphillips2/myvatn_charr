#=========================================================================================
#========== Preliminaries
#=========================================================================================

# load packages
library(tidyverse)
library(cowplot)
library(lemon)
library(nlme)
library(AICcmodavg)
library(demogR)

# import data 
site_data <- read_csv("model/site_data.csv") %>%
  mutate(stage = factor(stage,
                        levels = c("first",
                                   "second",
                                   "third",
                                   "adult"),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+")))


# import data_list
data_list <- read_rds("model/model_fits/primary/data_list.rds")

# import expanded model including separate surv sd's
fit_expanded <- read_rds("model/model_fits/expanded/fit.rds") 
fit_sum_expanded  <- read_csv("model/model_fits/expanded/fit_sum.csv")  

# import primary model with fixed surv sd's 
fit <- read_rds("model/model_fits/primary/fit.rds")                   
fit_sum <- read_csv("model/model_fits/primary/fit_sum.csv")

# reduced model with fixed sd's
fit_reduced <- read_rds("model/model_fits/reduced/fit.rds")           
fit_sum_reduced <- read_csv("model/model_fits/reduced/fit_sum.csv")

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(t = 1,
                                       r = 1,
                                       b = 1,
                                       l = 1),
                  legend.margin = margin(t = 0,
                                         r = 0,
                                         b = 0,
                                         l = -4),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black",family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,5,0,0)),
                  axis.title.x = element_text(margin = margin(5,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))

# year breaks
year_breaks <- c(1990, 2005, 2020)
year_limits <- c(1986, 2020.2)

# stage colors
stage_colors <- c("firebrick","dodgerblue","magenta4","goldenrod")
names(stage_colors) <- c("age 1","age 2","age 3","age 4+")

#=========================================================================================





#=========================================================================================
#========== Model comparison
#=========================================================================================

# posterior log-likelihood
fit_expanded_ll <- rstan::extract(fit_expanded, pars = "log_lik_sum") %>%
  as_tibble() %>%
  summarize(lo = quantile(log_lik_sum, probs = 0.16),
            mi = median(log_lik_sum),
            hi = quantile(log_lik_sum, probs = 0.85))
fit_ll <- rstan::extract(fit, pars = "log_lik_sum") %>%
  as_tibble() %>%
  summarize(lo = quantile(log_lik_sum, probs = 0.16),
            mi = median(log_lik_sum),
            hi = quantile(log_lik_sum, probs = 0.85))
fit_reduced_ll <- rstan::extract(fit_reduced, pars = "log_lik_sum") %>%
  as_tibble() %>%
  summarize(lo = quantile(log_lik_sum, probs = 0.16),
            mi = median(log_lik_sum),
            hi = quantile(log_lik_sum, probs = 0.85))

# examine surv sd's
fit_sum_expanded %>% filter(str_detect(.$var, "sls"))
fit_sum %>% filter(str_detect(.$var, "sls"))

#=========================================================================================





#=========================================================================================
#========== Population estimates: compare with catch data
#=========================================================================================

### Full model

# extract scaling parameter
k <- data_list$k

# extract detection probabilities
detect_prob_pars <- {fit_sum %>%
    filter(str_detect(.$var, "p\\["))}$var
detect_prob <- rstan::extract(fit, pars = detect_prob_pars) %>%
  parallel::mclapply(as_tibble) %>%
  bind_cols() %>%
  set_names(detect_prob_pars) %>%
  mutate(step = row_number()) %>%
  gather(var, val, -step) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+")))

# extract theta
theta <- rstan::extract(fit, pars = "theta") %>%
  parallel::mclapply(as_tibble) %>%
  bind_cols() %>%
  mutate(step = row_number())

# extract population density from MCMC
x_pars <- {fit_sum %>%
    filter(str_detect(.$var, "x\\["))}$var
x_full <- rstan::extract(fit, pars = x_pars) %>%
  parallel::mclapply(as_tibble) %>%
  bind_cols() %>%
  set_names(x_pars) %>%
  mutate(step = row_number()) %>%
  gather(var, val, -step) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         year = sort(unique(site_data$year))[time],
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+"))) 

# simulate prediction interval and summarize 90%
set.seed(3e2)
x_pred <- x_full %>%
  full_join(detect_prob %>%
              select(step, stage, val) %>%
              rename(p = val)) %>%
  full_join(theta %>%
              rename(theta = value)) %>%
  mutate(y_sim = rbinom(n = length(val), size =  1, prob = 1 - theta) * 
           rpois(n = length(val), k * p * val)) %>%
  group_by(stage, year) %>%
  summarize(lo = quantile(y_sim, probs = c(0.05)),
            mi = quantile(y_sim, probs = c(0.5)),
            hi = quantile(y_sim, probs = c(0.95))) %>%
  ungroup()

# create jitter manually (to be constent across panels)
set.seed(2e2)
site_data_jitter<- site_data %>%
  mutate(year_jitter = year + rnorm(n = length(x_pred$year), mean  = 0, sd = 0.2))


# plot annotation
labs <- x_full %>%
  tidyr::expand( stage) %>%
  mutate(x = mean(year_limits),
         y = 130)

# plot
p_catch_a <- ggplot(data = x_pred,
                  aes(x = year,
                      y = mi,
                      color = stage,
                      fill = stage))+
  facet_rep_wrap(~stage)+
  geom_text(data = labs,
            aes(label = stage,
                x = x,
                y = y),
            color = "black",
            size = 3.5,
            inherit.aes = F)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.1,
              linetype = 0)+
  geom_point(data = site_data_jitter,
             aes(x = year_jitter,
                 y = count),
             shape = 16,
             size = 0.5,
             alpha = 0.3)+
  geom_line(size = 0.6)+
  geom_line(data = site_data_jitter %>% 
              group_by(year, stage) %>%
              summarize(count = mean(count)),
            aes(x = year,
                y = count),
            linetype = "33",
            size = 0.5,
            color = "gray30")+
  scale_y_continuous("Survey catch",
                     trans = "log1p",
                     breaks = c(0, 5, 25, 125))+
  scale_x_continuous("",
                     breaks = year_breaks,
                     labels = NULL,
                     limits = year_limits + c(-0.5, 0.5))+
  scale_color_manual(values = stage_colors,
                     guide = "none")+
  scale_fill_manual(values = stage_colors,
                    guide = "none")+
  theme(plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(-0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
p_catch_a




### Reduced model
# extract scaling parameter
k <- data_list$k

# extract detection probabilities
detect_prob_pars_reduced <- {fit_sum_reduced %>%
    filter(str_detect(.$var, "p\\["))}$var
detect_prob_reduced <- rstan::extract(fit, pars = detect_prob_pars_reduced) %>%
  parallel::mclapply(as_tibble) %>%
  bind_cols() %>%
  set_names(detect_prob_pars_reduced) %>%
  mutate(step = row_number()) %>%
  gather(var, val, -step) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+")))

# extract theta
theta_reduced <- rstan::extract(fit_reduced, pars = "theta") %>%
  parallel::mclapply(as_tibble) %>%
  bind_cols() %>%
  mutate(step = row_number())

# extract population density from MCMC
x_pars_reduced <- {fit_sum_reduced %>%
    filter(str_detect(.$var, "x\\["))}$var
x_full_reduced <- rstan::extract(fit_reduced, pars = x_pars_reduced) %>%
  parallel::mclapply(as_tibble) %>%
  bind_cols() %>%
  set_names(x_pars_reduced) %>%
  mutate(step = row_number()) %>%
  gather(var, val, -step) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         year = sort(unique(site_data$year))[time],
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+"))) 

# simulate prediction interval and summarize 90%
set.seed(3e2)
x_pred_reduced <- x_full_reduced %>%
  full_join(detect_prob_reduced %>%
              select(step, stage, val) %>%
              rename(p = val)) %>%
  full_join(theta_reduced %>%
              rename(theta = value)) %>%
  mutate(y_sim = rbinom(n = length(val), size =  1, prob = 1 - theta) * 
           rpois(n = length(val), k * p * val)) %>%
  group_by(stage, year) %>%
  summarize(lo = quantile(y_sim, probs = c(0.05)),
            mi = quantile(y_sim, probs = c(0.5)),
            hi = quantile(y_sim, probs = c(0.95))) %>%
  ungroup()

# plot annotation
labs <- x_full_reduced %>%
  tidyr::expand( stage) %>%
  mutate(x = mean(year_limits),
         y = 130)

# plot
p_catch_b <- ggplot(data = x_pred_reduced,
                          aes(x = year,
                              y = mi,
                              color = stage,
                              fill = stage))+
  facet_rep_wrap(~stage)+
  geom_text(data = labs,
            aes(label = stage,
                x = x,
                y = y),
            color = "black",
            size = 3.5,
            inherit.aes = F)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.1,
              linetype = 0)+
  geom_point(data = site_data_jitter,
              aes(x = year_jitter,
                  y = count),
              shape = 16,
              size = 0.5,
              alpha = 0.3)+
  geom_line(size = 0.6)+
  geom_line(data = site_data_jitter %>% 
               group_by(year, stage) %>%
               summarize(count = mean(count)),
             aes(x = year,
                 y = count),
             linetype = "33",
             size = 0.5,
             color = "gray30")+
  scale_y_continuous("Survey catch",
                     trans = "log1p",
                     breaks = c(0, 5, 25, 125))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits+c(-0.5, 0.5))+
  scale_color_manual(values = stage_colors,
                     guide = "none")+
  scale_fill_manual(values = stage_colors,
                    guide = "none")+
  theme(plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(-0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
p_catch_b


# combine
p_catch <- plot_grid(NULL, p_catch_a, NULL, p_catch_b,
                    nrow = 4,
                    rel_heights = c(0.12, 1, 0.05, 1),
                    align = "h",
                    labels = c("",
                               "a: time-varying rates",
                               "",
                               "b: fixed rates"),
                    label_size = 10,
                    label_fontface = "plain",
                    hjust = c(0, 0, 0, 0),
                    vjust = c(0, -1.4, 0, -1.4))

# examine plot
p_catch

# export
# cairo_pdf(file = "analysis/figures/p_catch.pdf",
#           width = 3.5, height = 6, family = "Arial")
# p_catch
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Population estimates: relative density
#=========================================================================================

# extract scaling parameter
k <- data_list$k

# extract density estimate
x_fit <- fit_sum %>%
  filter(str_detect(.$var, "x\\[")) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         year = sort(unique(site_data$year))[time],
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+"))) %>%
  # scale by k
  mutate(lo = k * lo,
         mi = k * mi,
         hi = k * hi) %>%
  select(year, stage, lo, mi, hi)


# plot annotation
labs <- x_fit %>%
  filter(year == 2020) %>%
  mutate(x = 2020,
         y = mi)

# plot
p_dens <- ggplot(data = x_fit,
                 aes(x = year,
                     y = mi,
                     color = stage,
                     fill = stage))+
  geom_text(data = labs,
            aes(label = stage,
                x = x,
                y = y,
                color = stage),
            size = 3.5,
            inherit.aes = F,
            nudge_y = c(0, 0, 0, -3),
            nudge_x = c(3, 3, 3, -11))+
  geom_line(size = 0.5)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2,
              linetype = 0)+
  scale_y_continuous(Estimated~density~(`#`~station^{-1}),
                     trans = "log",
                     breaks = c(0.5, 5, 50, 500, 5000),
                     labels = c("0.5", "5", "50", "500", "5000"))+
  scale_x_continuous("Year",
                     breaks = c(1990, 2000, 2010, 2020),
                     limits = year_limits + c(0, 4))+
  scale_color_manual(values = stage_colors,
                     guide = "none")+
  scale_fill_manual(values = stage_colors,
                    guide = "none")+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))

# examine plot
p_dens

# export
# cairo_pdf(file = "analysis/figures/p_dens.pdf",
#           width = 3.5, height = 3, family = "Arial")
# p_dens
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Survival probability
#=========================================================================================

# logit survival
ls_fit <- fit_sum %>%
  filter(str_detect(.$var, "ls\\["),
         !str_detect(.$var, "sls\\[")) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         year = sort(unique(site_data$year))[time],
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+"))) 

# function to fit AR models for trends
ar_fit_fn <- function(stage_) {
  
  # extract relevant stage
  data_ = ls_fit %>% filter(stage == stage_)
  
  # define values for z-scoring
  mu_ = mean(data_$mi)
  sd_ = sd(data_$mi)
  
  # fit model
  m_ = gls(z ~ time,
           corAR1(form = ~ year),
           data  = data_ %>% 
             mutate(z = (mi - mu_) / sd_,
                    time = (year - mean(year)) / sd(year)))
  
  # summarize
  sum_ = summary(m_)
  
  # new data for prediction
  nd_ = data.frame(stage = stage_,
                   year = c(year = seq(min(data_$year), 
                                       max(data_$year),
                                       length.out = 100))) %>% 
    mutate(time = (year - mean(year)) / sd(year))
  
  # fitted values
  fit_ = predictSE.gls(m_, newdata = nd_, print.matrix = T)
  fit_ = cbind(nd_, fit_) %>%
    as_tibble() %>%
    mutate(mu = mu_,
           sd = sd_)
  
  # return
  return(list(model = m_,
              summary = sum_,
              fit = fit_))
}

# apply function to stages
ar_fit <- lapply(c("age 1","age 2","age 3","age 4+"),
                 ar_fit_fn) %>% 
  set_names(c("age 1","age 2","age 3","age 4+"))

# summary
lapply(ar_fit, function(x_){x_$summary})

# fitted values
ls_pred <- lapply(ar_fit, function(x_){x_$fit}) %>% 
  bind_rows() %>%
  mutate(fit = sd * fit + mu,
         se.fit = sd * fit + mu,
         stage = factor(stage,
                        levels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+")))

# plot annotation
labs <- x_full %>%
  tidyr::expand(stage) %>%
  mutate(x = mean(year_limits),
         y = 7)

surv_breaks <- c(0.01, 0.5, 0.99)
# plot
p_surv <- ggplot(data = ls_fit,
                 aes(x = year, 
                     y = mi,
                     color = stage))+
  facet_rep_wrap(~stage)+
  geom_hline(yintercept = 0, 
             linetype = 2, 
             size = 0.5,
             color = "gray50")+
  geom_text(data = labs,
            aes(label = stage,
                x = x,
                y = y),
            color = "black",
            size = 3.5,
            inherit.aes = F)+
  geom_ribbon(aes(ymin = lo,
                    ymax = hi,
                  fill = stage),
                linetype = 0,
                alpha = 0.2)+
  geom_line(size = 0.5)+
  geom_line(data = ls_pred,
            inherit.aes = F,
            aes(x = year,
                y = fit),
            size = 0.5)+
  scale_y_continuous("Survival probability",
                     breaks =log(surv_breaks / (1 - surv_breaks)),
                     limits = c(-7.2, 7.2),
                     labels = surv_breaks)+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits)+
  scale_color_manual(values = stage_colors,
                     guide = "none")+
  scale_fill_manual(values = stage_colors,
                    guide = "none")+
  theme(plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))

# examine plot
p_surv

# export
# cairo_pdf(file = "analysis/figures/p_surv.pdf",
#           width = 3.5, height = 3, family = "Arial")
# p_surv
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Recruitment
#=========================================================================================

# log recruitment
lr_fit <- fit_sum %>%
  filter(str_detect(.$var, "lr\\["),
         !str_detect(.$var, "slr\\[")) %>%
  mutate(time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         year = sort(unique(site_data$year))[time]) 


# define values for z-scoring
mu_rec = mean(lr_fit$mi)
sd_rec = sd(lr_fit$mi)

# fit model
rec_ar = gls(z ~ time,
             correlation = corAR1(form = ~year),
             data  = lr_fit %>% 
               mutate(z = (mi - mu_rec) / sd_rec,
                      time = (year - mean(year)) / sd(year)))

# summarize
summary(rec_ar)

# new data for prediction
rec_nd = data.frame(year = c(year = seq(min(lr_fit$year), 
                                     max(lr_fit$year),
                                     length.out = 100))) %>% 
  mutate(time = (year - mean(year)) / sd(year))

# fitted values
rec_pred = cbind(rec_nd,
                 predictSE.gls(rec_ar, newdata = rec_nd, print.matrix = T)) %>%
  as_tibble() %>%
  mutate(fit = sd_rec * fit + mu_rec,
         se.fit = sd_rec * se.fit + mu_rec)

# plot
p_rec <- ggplot(data = lr_fit,
                aes(x = year, 
                    y = mi))+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
                linetype = 0,
                alpha = 0.2,
                fill = "mediumblue")+
  geom_line(size = 0.5,
            color = "mediumblue")+
  geom_line(data = rec_pred,
            inherit.aes = F,
            aes(x = year,
                y = fit),
            size = 0.5)+
  scale_y_continuous(Recruitment~capita^{-1},
                     limits = c(0.5, 8.5),
                     breaks = log(c(3, 30, 300, 3000)),
                     labels = c(3, 30, 300, 3000))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))

# examine plot
p_rec

# export
# cairo_pdf(file = "analysis/figures/p_rec.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p_rec
# dev.off()


#=========================================================================================





#=========================================================================================
#========== Asymptotic growth rate
#=========================================================================================

# extract population density from MCMC (subset 2000)
# set.seed(1e3)
# dem_pars <- {fit_sum %>%
#     filter(str_detect(.$var, "s\\[") | str_detect(.$var, "r\\["),
#            !str_detect(.$var, "ls\\["),
#            !str_detect(.$var, "zs\\["),
#            !str_detect(.$var, "lr\\["),
#            !str_detect(.$var, "zr\\["))}$var
# dem_full <- rstan::extract(fit, pars = dem_pars) %>%
#   parallel::mclapply(as_tibble) %>%
#   bind_cols() %>%
#   set_names(dem_pars) %>%
#   sample_n(2000) %>%
#   mutate(step = row_number()) %>%
#   gather(var, val, -step) %>%
#   mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
#          age = ifelse(name == "r",
#                       1,
#                       strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))),
#          time = ifelse(name == "r",
#                        strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
#                        strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3]))),
#          year = sort(unique(site_data$year))[time])
# 
# # define matrix of zeros for storing values
# mat0 <- matrix(0, nrow = 4, ncol = 4)
# 
# # fill matrix and calculate asymptotic growth rate
# lambda_full <- dem_full %>%
#   split(.$step) %>%
#   parallel::mclapply(function(x_){
#     x_ %>% split(.$year) %>%
#       lapply(function(xx_){
#         v_ = xx_$val
#         mat_ <- mat0
#         mat_[1, 4] <- v_[1]
#         mat_[2, 1] <- v_[2]
#         mat_[3, 2] <- v_[3]
#         mat_[4, 3] <- v_[4]
#         mat_[4, 4] <- v_[5]
#         lam_ = eigen.analysis(mat_)$lambda1
#         return(tibble(year = unique(xx_$year),
#                       lam = lam_))
#       }) %>%
#       bind_rows()
#   }) %>% bind_rows()
# 
# write_csv(lambda_full, "analysis/lambda_full.csv")

lambda_full <- read_csv("analysis/lambda_full.csv")

# mean lambda
lambda_full %>%
  group_by(year) %>%
  mutate(step = row_number()) %>%
  group_by(step) %>%
  summarize(lam = prod(lam) ^ (1 / length(lam))) %>%
  ungroup() %>%
  summarize(lo = quantile(lam, probs = c(0.16)),
            mi = quantile(lam, probs = c(0.5)),
            hi = quantile(lam, probs = c(0.84))) %>%
  as.data.frame()

# summarize
lambda <- lambda_full %>%
  group_by(year) %>%
  summarize(lo = quantile(lam, probs = c(0.16)),
            mi = quantile(lam, probs = c(0.5)),
            hi = quantile(lam, probs = c(0.84)))


# define values for z-scoring
mu_lam = mean(log(lambda$mi))
sd_lam = sd(log(lambda$mi))

# fit model
lam_ar = gls(z ~ time,
             corAR1(form = ~year),
             data  = lambda %>% 
               mutate(z = (log(mi) - mu_lam) / sd_lam,
                      time = (year - mean(year)) / sd(year)))

# summarize
summary(lam_ar)

# new data forr prediction
lam_nd = data.frame(year = seq(min(lambda$year), 
                               max(lambda$year),
                               length.out = 100)) %>% 
  mutate(time = (year - mean(year)) / sd(year))

# fitted values
lam_fit = cbind(lam_nd,
                predictSE.gls(lam_ar, newdata = lam_nd, print.matrix = T)) %>%
  as_tibble() %>%
  mutate(fit = exp(sd_lam * fit + mu_lam),
         se.fit = exp(sd_lam * se.fit + mu_lam))


# plot
p_lam <- ggplot(data = lambda,
                aes(x = year, 
                    y = mi))+
  geom_hline(yintercept = 1, 
             linetype = 2, 
             size = 0.5,
             color = "gray50")+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              linetype = 0,
              alpha = 0.2,
              fill = "mediumblue")+
  geom_line(size = 0.5,
            color = "mediumblue")+
  geom_line(data = lam_fit,
            inherit.aes = F,
            aes(x = year,
                y = fit),
            size = 0.5)+
  scale_y_continuous(Population~growth~rate~(),
                     trans = "log",
                     breaks = c(1/3, 1, 3),
                     labels = c("1/3","1","3"),
                     limits = c(1/4, 3))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits)+
  theme(plot.margin = margin(t = 1,
                             r = 5,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
 
# examine plot
p_lam

# export
# cairo_pdf(file = "analysis/figures/p_lam.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p_lam
# dev.off()

#=========================================================================================




#=========================================================================================
#========== Ricker model
#=========================================================================================

# prepare data
ricker_data <- x_fit %>%
  group_by(year) %>%
  summarize(pop = sum(mi)) %>%
  full_join(lr_fit %>%
              select(year, mi) %>%
              rename(lpc_rec = mi)) %>%
  ungroup() %>%
  full_join(x_fit %>%
              filter(stage == "age 4+") %>%
              select(year, mi) %>%
              rename(adults = mi)) %>%
  mutate(rec = exp(lpc_rec) * adults) %>%
  na.omit()

# fit Ricker stock-recruitment curve on log scale
ricker_model <- nls(log(rec) ~ a + log(pop) - pop * (1/b),
         start = list(a = 1, b = 1),
         data = ricker_data)
summary(ricker_model)

# generated fitted valuees
ricker_fit <- tibble(pop = seq(min(ricker_data$pop),
                                   max(ricker_data$pop), 
                                   length.out = 1000))
ricker_fit$fit <- exp(predict(ricker_model, newdata = ricker_fit))

# plot
p_ricker <- ggplot(data = ricker_data,
                   aes(x = pop, 
                       y=rec))+
  geom_vline(xintercept = coef(ricker_model)["b"],
             linetype = 2,
             size = 0.5,
             color = "gray50")+
  geom_point(size = 1.75,
             alpha = 0.5)+
  geom_line(data = ricker_fit,
            size = 0.5,
            aes(x = pop,
                y = fit))+
  scale_y_continuous(Recruitment~(`#`~station^{-1}),
                     trans = "log",
                     breaks = c(30, 300, 3000))+
  scale_x_continuous(Population~density~(`#`~station^{-1}),
                     breaks = c(0, 2000, 4000, 6000))+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))


# examine plot
p_ricker

# export
# cairo_pdf(file = "analysis/figures/p_ricker.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p_ricker
# dev.off()

#=========================================================================================
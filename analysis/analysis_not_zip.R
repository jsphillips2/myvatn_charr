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
data <- read_csv("data/myvatn_char_clean_1986_2020.csv")
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


# import model fit
data_list <- read_rds(paste0("model/output/full/data_list.rds"))
fit <- read_rds(paste0("model/output/full/fit.rds"))
fit_sum <- read_csv(paste0("model/output/full/fit_sum.csv"))
fit_reduced <- read_rds(paste0("model/output/fixed_all/fit.rds"))
fit_sum_reduced <- read_csv(paste0("model/output/fixed_all/fit_sum.csv"))

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
year_limits <- c(1986, 2020)

# stage colors
stage_colors <- c("firebrick","dodgerblue","magenta4","goldenrod")
names(stage_colors) <- c("age 1","age 2","age 3","age 4+")

#=========================================================================================





#=========================================================================================
#========== Model comparison: log-likelihood
#=========================================================================================

# posterior log-likelihood
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

#=========================================================================================





#=========================================================================================
#========== Population estimates: compare with catch data
#=========================================================================================

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

# extract population density from MCMC (subset 2000)
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
         year = sort(unique(data$year))[time],
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+"))) 

# simulate prediction interval and summarize 90%
x_pred <- x_full %>%
  full_join(detect_prob %>%
              select(step, stage, val) %>%
              rename(p = val)) %>%
  mutate(y_sim = rpois(n = length(val), k * p * val)) %>%
  group_by(stage, year) %>%
  summarize(lo = quantile(y_sim, probs = c(0.05)),
            mi = quantile(y_sim, probs = c(0.5)),
            hi = quantile(y_sim, probs = c(0.95))) %>%
  ungroup()

# plot annotation
labs <- x_full %>%
  tidyr::expand( stage) %>%
  mutate(x = mean(year_limits),
         y = 130)

# plot
p_catch <- ggplot(data = x_pred,
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
              alpha = 0.2,
              linetype = 0)+
  geom_jitter(data = site_data,
              aes(x = year,
                  y = count),
              shape = 16,
              size = 0.5,
              alpha = 0.5)+
  geom_line(size = 0.5)+
  scale_y_continuous("Survey catch",
                     trans = "log1p",
                     breaks = c(0, 10, 100))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits)+
  scale_color_manual(values = stage_colors,
                     guide = F)+
  scale_fill_manual(values = stage_colors,
                    guide = F)+
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
p_catch

# cairo_pdf(file = "analysis/figures/p_catch_not_zip.pdf",
#           width = 3.5, height = 3, family = "Arial")
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
         year = sort(unique(data$year))[time],
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
                     breaks = c(0.5,  50,  5000),
                     labels = c("0.5",  "50",  "5000"))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits + c(0, 4))+
  scale_color_manual(values = stage_colors,
                     guide = F)+
  scale_fill_manual(values = stage_colors,
                    guide = F)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
p_dens

# cairo_pdf(file = "analysis/figures/p_dens_not_zip.pdf",
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
         year = sort(unique(data$year))[time],
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

# plot
p_surv <- ggplot(data = ls_fit,
                 aes(x = year, 
                     y = mi,
                     color = stage))+
  facet_rep_wrap(~stage)+
  geom_hline(yintercept = 0, 
             linetype = 2, 
             size = 0.3,
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
  scale_y_continuous("logit (survival probability)",
                     breaks = c(-5, 0, 5),
                     limits = c(-7.2, 8))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits)+
  scale_color_manual(values = stage_colors,
                     guide = F)+
  scale_fill_manual(values = stage_colors,
                    guide = F)+
  theme(plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
p_surv

# cairo_pdf(file = "analysis/figures/p_surv_not_zip.pdf",
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
         year = sort(unique(data$year))[time]) 


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
  scale_y_continuous(log~(recruitment~capita^{-1}),
                     limits = c(-1, 9),
                     breaks = c(-1, 2, 5, 8))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
p_rec

# cairo_pdf(file = "analysis/figures/p_rec_not_zip.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p_rec
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Asymptotic growth rate
#=========================================================================================

# # extract population density from MCMC (subset 2000)
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
#          year = sort(unique(data$year))[time])
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
# write_csv(lambda_full, "analysis/lambda_full_not_zip.csv")

# lambda_full <- read_csv("analysis/lambda_full_not_zip.csv")

# mean lambda
lambda_full %>%
  group_by(year) %>%
  mutate(step = row_number()) %>%
  group_by(step) %>%
  summarize(lam = prod(lam) ^ (1 / length(lam))) %>%
  ungroup() %>%
  summarize(lo = quantile(lam, probs = c(0.16)),
            mi = quantile(lam, probs = c(0.5)),
            hi = quantile(lam, probs = c(0.84)))

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
             size = 0.3,
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
  scale_y_continuous("Population growth rate",
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
  
p_lam

# cairo_pdf(file = "analysis/figures/p_lam_not_zip.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# p_lam
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Reduced model population estimates: compare with catch data
#=========================================================================================

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


# extract population density from MCMC (subset 2000)
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
         year = sort(unique(data$year))[time],
         stage = factor(age,
                        levels = c(1,2,3,4),
                        labels = c("age 1",
                                   "age 2",
                                   "age 3",
                                   "age 4+"))) 

# simulate prediction interval and summarize 90%
x_pred_reduced <- x_full_reduced %>%
  full_join(detect_prob_reduced %>%
              select(step, stage, val) %>%
              rename(p = val)) %>%
  mutate(y_sim = rpois(n = length(val), k * p * val)) %>%
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
p_catch_reduced <- ggplot(data = x_pred_reduced,
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
              alpha = 0.2,
              linetype = 0)+
  geom_jitter(data = site_data,
              aes(x = year,
                  y = count),
              shape = 16,
              size = 0.5,
              alpha = 0.5)+
  geom_line(size = 0.5)+
  scale_y_continuous("Survey catch",
                     trans = "log1p",
                     breaks = c(0, 10, 100))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits)+
  scale_color_manual(values = stage_colors,
                     guide = F)+
  scale_fill_manual(values = stage_colors,
                    guide = F)+
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
p_catch_reduced

# cairo_pdf(file = "analysis/figures/p_catch_reduced_not_zip.pdf",
#           width = 3.5, height = 3, family = "Arial")
# p_catch_reduced
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Reduced population estimates: relative density
#=========================================================================================

# extract scaling parameter
k <- data_list$k

# extract density estimate
x_fit_reduced <- fit_sum_reduced %>%
  filter(str_detect(.$var, "x\\[")) %>%
  mutate(age = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         year = sort(unique(data$year))[time],
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
labs <- x_fit_reduced %>%
  filter(year == 2020) %>%
  mutate(x = 2020,
         y = mi)

# plot
p_dens_reduced <- ggplot(data = x_fit_reduced,
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
            nudge_y = c(0, 0, 0, 0),
            nudge_x = c(3.1, 3.1, 3.1, 3.6))+
  geom_line(size = 0.5)+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2,
              linetype = 0)+
  scale_y_continuous(Estimated~density~(`#`~station^{-1}),
                     trans = "log",
                     breaks = c(5, 25, 125),
                     limits = c(2,200))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits + c(0, 5))+
  scale_color_manual(values = stage_colors,
                     guide = F)+
  scale_fill_manual(values = stage_colors,
                    guide = F)+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))
p_dens_reduced

# cairo_pdf(file = "analysis/figures/p_dens_reduced_not_zip.pdf",
#           width = 3.5, height = 3, family = "Arial")
# p_dens_reduced
# dev.off()

#=========================================================================================

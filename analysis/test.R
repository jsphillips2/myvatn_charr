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
data <- read_csv("data/myvatn_char_clean.csv")
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
fit <- read_rds(paste0("model/output/full_zip/fit.rds"))
fit_sum <- read_csv(paste0("model/output/full_zip/fit_sum.csv"))

# set theme
theme_set(theme_bw() %+replace%
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  plot.margin = margin(1,1,1,1),
                  legend.margin = margin(0,0,0,-4),
                  legend.text = element_text(size = 8),
                  axis.text = element_text(size = 10, color="black",family = "sans"),
                  axis.title = element_text(size =10),
                  axis.title.y = element_text(angle = 90, margin=margin(0,5,0,0)),
                  axis.title.x = element_text(margin = margin(5,0,0,0)),
                  panel.spacing = unit(0.1, "lines"),
                  axis.ticks = element_line(size = 0.25)))

# year breaks
year_breaks <- c(1985, 2000, 2015)
year_limits <- c(1985, 2017)

# stage colors
stage_colors <- c("firebrick","dodgerblue","magenta4","goldenrod")
names(stage_colors) <- c("age 1","age 2","age 3","age 4+")

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

# extract theta
theta <- rstan::extract(fit, pars = "theta") %>%
  parallel::mclapply(as_tibble) %>%
  bind_cols() %>%
  mutate(step = row_number())

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
  full_join(theta %>%
              rename(theta = value)) %>%
  mutate(y_sim = rbinom(n = length(val), size =  1, prob = 1 - theta) * 
                    rpois(n = length(val), k * p * val)) %>%
  group_by(stage, year) %>%
  summarize(lo = quantile(y_sim, probs = c(0.05)),
            mi = quantile(y_sim, probs = c(0.5)),
            hi = quantile(y_sim, probs = c(0.95))) %>%
  ungroup()

# plot annotation
labs <- x_full %>%
  tidyr::expand(stage) %>%
  mutate(x = mean(year_limits) + 1.5,
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
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))+
  coord_capped_cart(left = "both", 
                    bottom="both")
p_catch
# cairo_pdf(file = "analysis/p_catch_zip.pdf",
#           width = 3.5, height = 3, family = "Arial")
# p_catch
# dev.off()
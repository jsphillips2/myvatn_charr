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
age_data <- read_csv("analysis/age_data.csv")


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
year_limits <- c(1986, 2020)

# stage colors
stage_colors <- c("firebrick","dodgerblue","magenta4","goldenrod")
names(stage_colors) <- c("age 1","age 2","age 3","age 4+")

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
labs <- x_fit_reduced %>%
  filter(year == 2020) %>%
  mutate(x = 2020,
         y = mi)

# plot
s_dens_reduced <- ggplot(data = x_fit_reduced,
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
                     limits = c(3.5,200))+
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

# examine
s_dens_reduced

# export
# cairo_pdf(file = "analysis/figures/s_dens_reduced.pdf",
#           width = 3.5, height = 3, family = "Arial")
# s_dens_reduced
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Catch by site
#=========================================================================================

# plot annotation
labs <- x_full %>%
  tidyr::expand(stage) %>%
  mutate(x = mean(year_limits) + 1.5,
         y = 130)

# plot
s_sites <- ggplot(data = site_data,
                  aes(x = year,
                      y = count,
                      color = factor(area + 1)))+
  facet_rep_wrap(~stage)+
  geom_text(data = labs,
            aes(label = stage,
                x = x,
                y = y),
            color = "black",
            size = 3.5,
            inherit.aes = F)+
  geom_line(size = 0.3,
            alpha = 0.75)+
  scale_y_continuous("Survey catch",
                     trans = "log1p",
                     breaks = c(0, 10, 100))+
  scale_x_continuous("Year",
                     breaks = year_breaks,
                     limits = year_limits)+
  scale_fill_manual(values = stage_colors,
                    guide = F)+
  scale_color_viridis_d("")+
  theme(legend.position = "top",
        plot.margin = margin(t = 1,
                             r = 10,
                             b = 1,
                             l = 1),
        panel.border = element_blank(),
        panel.spacing.x = unit(-0.5, "lines"),
        panel.spacing.y = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))

# examine
s_sites

# export
# cairo_pdf(file = "analysis/figures/s_sites.pdf",
#           width = 3.5, height = 4, family = "Arial")
# s_sites
# dev.off()

#=========================================================================================





#=========================================================================================
#========== Catch by cohort
#=========================================================================================

# prepare cohort data
# fill in missing ages with estimated age and drop June data
cohorts <- age_data %>%
  mutate(cohort = year - est_age,
         cohort = cohort - min(cohort) + 1) %>%
  group_by(cohort) %>%
  {full_join(., expand(., cohort = min(cohort):max(cohort), est_age = min(est_age):max(est_age)))} %>%
  mutate(
    count = ifelse(is.na(count)==T | count == 0, 1L, count) 
  ) %>%
  ungroup()

# plot
s_cohort <- ggplot(data = cohorts,
                   aes(x = est_age,
                       y = count,
                       group = cohort))+
  geom_line(size = 0.3,
            alpha = 0.75)+
  scale_y_continuous("Catch by cohort",
                     trans="log",
                     breaks=c(1,10,100, 1000),
                     limits = c(1, 1000))+
  scale_x_continuous("Cohort age",
                     breaks = c(1, 3, 5, 7, 9, 11))+
  theme(panel.border = element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_blank(),
        axis.line.x = element_line(size = 0.25),
        axis.line.y = element_line(size = 0.25))

# examine
s_cohort  

# export
# cairo_pdf(file = "analysis/figures/s_cohort.pdf",
#           width = 3.5, height = 2.5, family = "Arial")
# s_cohort
# dev.off()
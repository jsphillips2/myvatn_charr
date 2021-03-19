#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# import data
data<-read_csv("data/myvatn_char_clean_1986_2020.csv")

# set theme
theme_set(theme_bw() %+replace% 
            theme(panel.grid = element_blank(),
                  strip.background = element_blank(),
                  legend.margin = margin(0,0,0,0),
                  strip.text = element_text(size=10),
                  legend.text = element_text(size=10),
                  axis.text=element_text(size=10, color="black"),
                  axis.title.y=element_text(angle = 90 ,margin=margin(0,15,0,0)),
                  axis.title.x=element_text(margin=margin(15,0,0,0))))





#==========
#========== Prepare Data
#==========

clean_data <-  data %>%
  #Fill gaps in "Est.age"
  mutate(est_age = ifelse(is.na(age)==F, age, est_age)) %>%
  #Remove June data & samples without ages
  filter(month != 6, is.na(est_age)==F) %>%
  #Define Stage
  mutate(stage = factor(ifelse(est_age == 1, "first", 
                               ifelse(est_age==2, "second", 
                                      ifelse(est_age==3, "third", "adult"))),
                        levels = c("first","second","third","adult")),
         area = ifelse(area_1 == 11, 7 ,
                       ifelse(area_1 == 12, 8, area_1)))

#count number of fish per stage per year
site_data = clean_data %>%
  group_by(year, stage, area) %>%
  #Count number of fish
  summarize(
    count = length(length) 
  ) %>%
  #Merge with stage_year to keep track of 0's
  full_join(clean_data %>% expand(stage, year, area)) %>% 
  #Repalce NA's with 1's
  #This only affects data for 1st year fish, which aren't used in the model
  mutate(
    count = ifelse(is.na(count)==T, 0L, count) 
  ) %>%
  ungroup()
# write_csv(site_data, "model/site_data.csv")

# plot
ggplot(data = site_data,
       aes(year, count))+
  facet_wrap(~stage)+
  geom_line(aes(group = area), size = 0.3, alpha = 0.5)+
  geom_line(data = site_data %>% group_by(stage, year) %>% 
              summarize(count = mean(count,na.rm=T)),
            size = 0.9)+
  scale_x_continuous("Year", limits=c(1986,2020))+
  scale_y_continuous("Abundance by site",trans="log1p",breaks=c(1,10,100))
  
# define indices 
n_stages <- length(unique(site_data$stage))
n_sites <- length(unique(site_data$area))
n_years <- length(unique(site_data$year))

# data array
y <- array({site_data %>% arrange(year, area, stage)}$count,
           dim = c(n_stages, 
                   n_sites,
                   n_years))

# arrange values for each stage x year combination in descending order
n_zero <- matrix(0, nrow = n_stages, ncol = n_years)
n_nonzero <- matrix(0, nrow = n_stages, ncol = n_years)
for (t in 1 : n_years) {
  for (i in 1 : n_stages) {
    y[i, , t] = sort(y[i, ,t], decreasing = T)
    n_zero[i, t] = length(y[i, ,t][y[i, ,t] == 0])
    n_nonzero[i, t] = length(y[i, ,t][y[i, ,t] > 0])
  }
}

# check
plot(colSums(y[1,,]) /n_sites , type = "l")
plot(colSums(y[2,,]) /n_sites , type = "l")
plot(colSums(y[3,,]) /n_sites , type = "l")
plot(colSums(y[4,,]) /n_sites , type = "l")

plot(count ~ year, 
     type = "l",
     data = 
       site_data %>% filter(stage == "first") %>% group_by(year) %>% summarize(count = sum(count)))
plot(count ~ year, 
     type = "l",
     data = 
       site_data %>% filter(stage == "second") %>% group_by(year) %>% summarize(count = sum(count)))
plot(count ~ year, 
     type = "l",
     data = 
       site_data %>% filter(stage == "third") %>% group_by(year) %>% summarize(count = sum(count)))
plot(count ~ year, 
     type = "l",
     data = 
       site_data %>% filter(stage == "adult") %>% group_by(year) %>% summarize(count = sum(count)))


# priors
xp <- matrix(c(20,       20,
               0,        5,
               0,        1,
               2,        2,
               1.5,      1,
               1.5,      1,
               2,        2),
             nrow = 7,
             byrow = T)

# scaling for avoid excessively large values for x during fitting
k <- 100

# time varying rates indicator
tvr <- c(0, 0)

# package data
data_list <- list(n_stages = n_stages,
                  n_sites = n_sites,
                  n_years = n_years,
                  y = y,
                  xp = xp,
                  tvr = tvr,
                  k = k,
                  n_zero = n_zero,
                  n_nonzero = n_nonzero)




#==========
#========== Fit model
#==========

# MCMC specifications (for testing)
# chains <- 1
# iter <- 100
# adapt_delta <- 0.8
# max_treedepth <- 13

# MCMC specifications
# chains <- 4
# iter <- 2000
# adapt_delta <- 0.95
# max_treedepth <- 14

# # fit model
# start_time <- Sys.time()
# fit <- stan(file = "model/charr_model_zip.stan",
#             data = data_list,
#             seed=2e3,
#             chains = chains,
#             iter = iter,
#             control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))
# end_time <- Sys.time()
# end_time - start_time
# 
# 
# 
# fit_sum <- rstan::summary(fit, probs = c(0.16,0.50,0.84))$summary %>%
#     as.data.frame() %>%
#     rownames_to_column() %>%
#     as_tibble() %>%
#     rename(var = rowname,
#            lo = `16%`,
#            mi = `50%`,
#            hi = `84%`) %>%
#     select(var, lo, mi, hi, n_eff, Rhat)
# 
# # export
# name <- if (tvr[1] == 1 && tvr[2] == 1) {"full"} else {
#           if (tvr[1] == 0 && tvr[2] == 1) {"fixed_r"} else {
#             if (tvr[1] == 1 && tvr[2] == 0) {"fixed_s"} else {"fixed_all"}}}
# write_rds(data_list, paste0("model/output/",name,"_zip/data_list.rds"))
# write_rds(fit, paste0("model/output/",name,"_zip/fit.rds"))
# write_csv(fit_sum, paste0("model/output/",name,"_zip/fit_sum.csv"))





#==========
#========== Check fit
#==========

# name <- if (tvr[1] == 1 && tvr[2] == 1) {"full"} else {
#           if (tvr[1] == 0 && tvr[2] == 1) {"fixed_r"} else {
#             if (tvr[1] == 1 && tvr[2] == 0) {"fixed_s"} else {"fixed_all"}}}
# fit <- read_rds(paste0("model/output/",name,"/fit.rds"))
# fit_sum <- read_csv(paste0("model/output/",name,"/fit_sum.csv"))

# fit <- read_rds(paste0("model/output/full_zip/fit.rds"))
# fit_sum <- read_csv(paste0("model/output/full_zip/fit_sum.csv"))

x_fit <- fit_sum %>%
  filter(str_detect(.$var, "x\\[")) %>%
  mutate(stage = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3]))) 

x_fit %>%
  ggplot(aes(time, mi))+
  facet_wrap(~stage, scales = "free_y")+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2)


s_fit <- fit_sum %>%
  filter(str_detect(.$var, "s\\["),
         !str_detect(.$var, "ls\\["),
         !str_detect(.$var, "zs\\[")) %>%
  mutate(stage = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3]))) 

s_fit %>%
  ggplot(aes(time, mi))+
  facet_wrap(~stage)+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2)


r_fit <- fit_sum %>%
  filter(str_detect(.$var, "r\\["),
         !str_detect(.$var, "lr\\["),
         !str_detect(.$var, "zr\\[")) %>%
  mutate(time = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) 

r_fit %>%
  ggplot(aes(time, mi))+
  geom_line()+
  geom_ribbon(aes(ymin = lo,
                  ymax = hi),
              alpha = 0.2)+
  scale_y_continuous(trans = "log")


ps <- fit_sum %>%
  filter(str_detect(.$var, "p\\[")) %>%
  mutate(stage = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         stage = factor(stage,
                        levels = c(1,2,3,4),
                        labels = c("first",
                                   "second",
                                   "third",
                                   "adult")))

test <- x_fit %>% 
  mutate(stage = factor(stage,
                        levels = c(1,2,3,4),
                        labels = c("first",
                                   "second",
                                   "third",
                                   "adult")),
         year = 1985 + time) %>%
  full_join(ps %>% select(stage, mi) %>% rename(p = mi)) %>%
  mutate(lo = k * lo * p,
         mi = k * mi * p,
         hi = k * hi * p)

ggplot(data = site_data ,
       aes(year, count))+
  facet_wrap(~stage)+
  geom_line(aes(group = area), size = 0.3, alpha = 0.5)+
  geom_line(data = test,
            aes(y = mi),
            size = 0.9)+
  geom_ribbon(data = test,
              aes(ymin = lo,
                  ymax =hi,
                  x = year),
              inherit.aes = F,
              alpha = 0.2)+
  scale_x_continuous("Year", limits=c(1986,2020))+
  scale_y_continuous(trans = "log1p")

fit_sum %>%
  filter(str_detect(.$var, "slr"))

fit_sum %>%
  filter(str_detect(.$var, "sls"))

fit_sum %>%
  filter(str_detect(.$var, "p"),!str_detect(.$var, "lp"))

fit_sum %>%
  filter(var == "theta")

fit_sum %>%
  filter(str_detect(.$var, "ls0"))

x_pars <- {fit_sum %>%
  filter(str_detect(.$var, "x\\["))}$var
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
site_data <-read_csv("model/site_data.csv") %>%
  mutate(stage = factor(stage,
                        levels = c("first","second","third","adult")))

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
tvr <- c(1, 1) # full 
# tvr <- c(0, 0) # both fixed

# slsm <- c(1, 2, 3, 4) # separate sd's
slsm <- c(1, 1, 1, 1) # single sd

# package data
data_list <- list(n_stages = n_stages,
                  n_sites = n_sites,
                  n_years = n_years,
                  y = y,
                  xp = xp,
                  tvr = tvr,
                  slsm = slsm,
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
# fit <- stan(file = "charr_model.stan",
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
# name <- {if (tvr[1] == 1 && tvr[2] == 1 && max(slsm) == 4) {"expanded"} else 
#           if (tvr[1] == 1 && tvr[2] == 1 && max(slsm) == 1) {"primary"} else "reduced"}
# write_rds(data_list, paste0("model/model_fits/",name, "/data_list.rds"))
# write_rds(fit, paste0("model/model_fits/", name,"/fit.rds"))
# write_csv(fit_sum, paste0("model/model_fits/", name,"/fit_sum.csv"))

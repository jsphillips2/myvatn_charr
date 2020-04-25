#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(runjags)

# import data
data = read_csv("data/myvatn_char_clean.csv")

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
#========== Clean and examine data
#==========

# trim data
# define stages and aera
data_clean <- data %>%
  #Fill gaps in "Est.age"
  mutate(est_age = ifelse(is.na(age)==F, age, est_age)
  ) %>%
  #Remove June data & samples without ages
  filter(month != 6, is.na(est_age)==F) %>%
  #Define Stage
  mutate(stage = factor(ifelse(est_age == 1, "first", 
                               ifelse(est_age==2, "second", 
                                      ifelse(est_age==3, "third", "adult"))),
                        levels = c("first","second","third","adult")),
         area = ifelse(area_1 == 11, 7 ,
                       ifelse(area_1 == 12, 8, area_1)))



# define coutns by site
site_counts <-  data_clean %>%
  group_by(year, stage, area) %>%
  #Count number of fish
  summarize(
    count = length(length) 
  ) %>%
  #Merge with stage_year to keep track of 0's
  full_join(data_clean %>% expand(stage, year, area)) %>% 
  #Repalce NA's with 1's
  #This only affects data for 1st year fish, which aren't used in the model
  mutate(
    count = ifelse(is.na(count)==T, 0L, count) 
  ) %>%
  ungroup()





#==========
#========== Prepare data
#==========

# prepare data
data_prep <- site_counts %>%
  filter(stage != "first") %>%
  arrange(year, area, stage) %>%
  mutate(year_n = as.numeric(factor(year)) - min(as.numeric(factor(year))) + 1,
         area_n = as.numeric(factor(area)) - min(as.numeric(factor(area))) + 1,
         stage_n = as.numeric(factor(stage, 
                                     levels = c("second",
                                                "third",
                                                "adult"))),
         nu_group = as.numeric(factor(paste(year_n, stage_n))),
         phi_group = as.numeric(factor(paste(stage_n, area)))) %>%
  arrange(year, area, stage)

# define groupings
lambda_nu <- {data_prep %>%
  group_by(nu_group) %>%
  summarize(stage_n = unique(stage_n))}$stage_n

phi_y <- data_prep$phi_group

nu_y <- data_prep$nu_group

# lambda priors
data_prep %>%
  group_by(year, stage) %>%
  summarize(count = sum(count)) %>%
  group_by(stage) %>%
  summarize(count = mean(count))
lp <- matrix(c(1.5, 1.5/5000),
             nrow = 1, ncol = 2,
             byrow = T) 

# phi priors
pp <- matrix(c(1.5, 1.5/0.5 - 1.5),
             nrow = 1, ncol = 2,
             byrow = T) 

# sigma priors
rl <- matrix(c(1.5, .5/2),
             nrow = 1, ncol = 2,
             byrow = T) 

# indeces
l <- max(lambda_nu)
p <- max(phi_y)
n <- max(nu_y)

# counts
y <- data_prep$count

# package data
data_list <- list(l = l,
                  p = p,
                  n = n,
                  rl = rl,
                  lambda_nu = lambda_nu,
                  phi_y = phi_y,
                  nu_y = nu_y,
                  y = y,
                  lp = lp,
                  pp = pp)

# set seed
seed = 2e4

nu_start <- data_prep %>%
  group_by(nu_group, stage) %>%
  summarize(count = max(count)) %>%
  left_join(data_prep %>%
              expand(stage) %>%
              mutate(p = pp[,1] / (pp[,1] + pp[,2]))) %>%
  mutate(count_max = (count + 1)/p)

# function for initial values
init_fn <- function(){
  list(lambda = runif(l, 
                      0.5 * (lp[,1] / lp[,2]), 
                      1.5 * (lp[,1] / lp[,2])),
       r = runif(1, 
                      0.5 * (rl[,1] / rl[,2]), 
                      1.5 * (rl[,1] / rl[,2])),
       phi = runif(p, 0, 1),
       nu = as.integer(runif(n,
                             nu_start$count,
                             nu_start$count_max)),
       .RNG.name="base::Wichmann-Hill", .RNG.seed=seed)
}





### Fit model

# file path for model
model_path <- "population_estimate/char_pop_est.txt"

# variables to monitor
monitor <- c("lambda","nu","phi","r","nu_l")

# MCMC specifications (for testing)
# n.chains <- 1
# adapt <- 1
# burnin <- 1
# sample <- 1

# MCMC specifications
# n.chains <- 4
# adapt <- 40000
# burnin <- 40000
# sample <- 2000

# fit
# fit <- run.jags(model = model_path, data = data_list, inits = init_fn,
#                 n.chains = n.chains, burnin = burnin, sample = sample,
#                 adapt = adapt, monitor = monitor, thin = 20)


# check diagnostic
psrf_check <- coda::gelman.diag(fit$mcmc)
psrf_clean <- psrf_check$psrf %>%
  as_tibble() %>%
  mutate(rowname =  rownames(psrf_check$psrf)) %>%
  rename(est  = `Point est.`,
         upp = `Upper C.I.`) %>%
  select(rowname, est, upp) %>%
  arrange(-est)


# summarize
options(mc.cores = parallel::detectCores()-2)
fit_sum <- fit$mcmc %>%
  parallel::mclapply(function(x_){
    d_ = x_ %>%
      as_tibble()
    return(d_)
  }) %>%
  bind_rows() %>%
  gather(rowname, val) %>%
  group_by(rowname) %>%
  summarize(mean = mean(val),
            sd = sd(val),
            lower68 = quantile(val, probs = 0.16),
            median = quantile(val, probs = 0.5),
            upper68 = quantile(val, probs = 0.84)) 


nu_fit <- fit_sum %>%
  filter(str_detect(.$rowname, "nu\\[")) %>%
  mutate(nu_group = strsplit(rowname, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(mean, sd, lower68, median, upper68, nu_group) %>%
  full_join(data_prep %>% 
              expand(nesting(year, stage, nu_group))) %>%
  arrange(year, stage)


nu_fit %>% 
  ggplot(aes(year, median, color = stage, fill = stage))+
  geom_line()+
  geom_ribbon(aes(ymin = lower68, ymax = upper68),
              linetype = 0, alpha = 0.2)+
  scale_color_manual(values=c("gray50","dodgerblue2","firebrick2","black"))+
  scale_fill_manual(values=c("gray50","dodgerblue2","firebrick2","black"))+
  scale_y_continuous("Char Count",trans="log")


fit_sum %>%
  filter(str_detect(.$rowname, "phi\\["))

fit_sum %>%
  filter(rowname == "lambda[1]")

fit_sum %>%
  filter(rowname == "lambda[2]")

fit_sum %>%
  filter(rowname == "lambda[3]")


bayesplot::mcmc_trace(x = fit$mcmc, pars = "lambda[1]")

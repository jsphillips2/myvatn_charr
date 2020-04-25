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
  mutate(area_n = as.numeric(factor(area)) - min(as.numeric(factor(area))),
         stage_n = as.numeric(factor(stage,
                                     levels = c("second",
                                                "third",
                                                "adult"))) - 1,
         year_n = as.numeric(factor(year)) - min(as.numeric(factor(year)))) %>%
  arrange(area, stage, year) %>%
  mutate(area_stage_n = area_n * stage_n,
         area_year_n = area_n * year_n,
         stage_year_n = stage_n * year_n,
         area_stage_year_n = area_n * stage_n * year_n,
         l_group = as.numeric(factor(paste(stage_n, year_n))),
         n_group = as.numeric(factor(paste(stage_n, year_n))),
         p_group = as.numeric(factor(paste(area_n,
                                           stage_n,
                                           area_stage_n
         ))))



  

# extract observations
y <- data_prep$count %>% round(0)



### Define model matrix for detection probabilites

# select variables
x_prep <- data_prep %>%
  select(area_n,
         stage_n,
         area_stage_n) %>%
  unique()

# calculate shifting of levels by column 
x_shift <- c(1, 1  + cumsum(x_prep %>% apply(2, function(x_){max(unique(x_))})))

# function for shifting levels by column (leaving 0's as 'intercepts')
shift_fn <- function(x_, shift_) {ifelse(x_ == 0, 1, shift_)}

# shift levels by column
x <- lapply(1:ncol(x_prep), function(n) {
  data.frame(v = x_prep[, n]  + shift_fn(x_prep[, n], x_shift[n]))
}) %>%
  bind_cols() %>%
  as.matrix()

# indeces for model matrix
n_j = ncol(x)
n_k = nrow(x)
n_beta <- max(x)

### Define mappings

# basin abundance (lambda) to station abundance (nu) 
lam_nu <- {data_prep %>% 
    group_by(n_group) %>%
    summarize(l_group = unique(l_group)) %>% 
    ungroup() %>%
    arrange(n_group)}$l_group
n_nu <- length(lam_nu)
n_lam <- max(lam_nu)

# station abundance (nu) to trap observations (y)
nu_y <- {data_prep %>% 
    select(n_group)}$n_group

# detection probability (phi) to trap observations (y)
phi_y <- {data_prep %>% 
    select(p_group)}$p_group



### Priors
u_lam_r <- c(1000, 500) # basin abundance mean
s_lam_r  <- c(500, 50) # basin abundance sd 
u_beta_r <- c(-5, rep(0, n_beta-1)) # coefficient mean (separate value for intercept) 
s_beta_r <- c(2, rep(0.1, n_beta-1)) # coefficient sd (separate value for intercept)



### Package data
data_list = list(y = y, # observed trap counts
                 x = x, # model matrix for detection probabilities 
                 n_j = n_j, # number of columns in model matrix
                 n_k = n_k, # number of rows in model matrix 
                 n_beta = n_beta, # number of coefficients 
                 n_nu = n_nu, # number of station estimates 
                 n_lam = n_lam, # number of basin estimates 
                 lam_nu = lam_nu, # mapping basin to station abundance 
                 nu_y = nu_y, # mapping station abundance to observed trap counts
                 phi_y = phi_y, # mapping of detection probability to observed trap counts
                 u_lam_r = u_lam_r, # basin abundance mean
                 s_lam_r = s_lam_r, # basin abundance sd
                 u_beta_r = u_beta_r, # coefficient mean
                 s_beta_r = s_beta_r # coefficient sd
)





#==========
#========== Fit model
#==========



### Initial values

# basin abundance (mean across stations)
lp <- {data_prep %>%
    group_by(l_group) %>%
    summarize(count = mean(count) + 1)}$count

# station abundance (max across traps)
nmax <- {data_prep %>%
    group_by(n_group) %>%
    summarize(nmax = max(count))}$nmax

# set seed
seed = 2e4

# function for initial values
init_fn <- function(){
  list(u_lam = runif(1, 
                     u_lam_r[1] - abs(u_lam_r[1]/2), 
                     u_lam_r[1] + abs(u_lam_r[1]/2)),
       s_lam = runif(1, 
                     s_lam_r[1] - abs(s_lam_r[1]/2), 
                     s_lam_r[1] + abs(s_lam_r[1]/2)),
       beta = runif(n_beta, u_beta_r - s_beta_r, u_beta_r + s_beta_r),
       lam = runif(lp, 
                   (lp - 0.5*lp)/(exp(-2)/(1+exp(-2))), 
                   (lp + 0.5*lp)/(exp(-2)/(1+exp(-2)))),
       nu = as.integer(round(runif(length(nmax), 
                                   nmax/(exp(-2)/(1+exp(-2))), 
                                   1.5*(nmax + 1)/(exp(-2)/(1+exp(-2)))))),
       .RNG.name="base::Wichmann-Hill", .RNG.seed=seed)
}

# file path for model
model_path <- "population_estimate/pop_est.txt"

# variables to monitor
monitor <- c("u_lam","s_lam","lam","nu","beta","phi")

# MCMC specifications (for testing)
# n.chains <- 1
# adapt <- 1
# burnin <- 1
# sample <- 1

# MCMC specifications
# n.chains <- 4
# adapt <- 4000
# burnin <- 4000
# sample <- 4000

# fit
fit <- run.jags(model = model_path, data = data_list, inits = init_fn,
                n.chains = n.chains, burnin = burnin, sample = sample,
                adapt = adapt, monitor = monitor)

# export fit
# saveRDS(fit, paste0("analysis/population_estimate/fit.rds"))

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





phi_fit <- fit_sum %>%
  filter(str_detect(.$rowname, "phi\\[")) %>%
  mutate(p_group = strsplit(rowname, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(mean, sd, lower68, median, upper68, p_group) %>%
  full_join(data_prep %>% 
              expand(nesting(year, stage, area, p_group))) %>%
  arrange(year, stage)




phi_fit %>%
  ggplot(aes(area, mean, group = area, color = stage))+
  geom_point()+
  scale_color_manual(values=c("dodgerblue2","firebrick2"," black"))


fit_sum %>%
  filter(str_detect(.$rowname, "lam\\["))

fit_sum %>%
  filter(str_detect(.$rowname, "u_lam"))

fit_sum %>%
  filter(str_detect(.$rowname, "s_lam"))


bayesplot::mcmc_trace(x = fit$mcmc, pars = "phi[3]")






lam_fit <- fit_sum %>%
  filter(str_detect(.$rowname, "lam\\[")) %>%
  mutate(l_group = strsplit(rowname, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(mean, sd, lower68, median, upper68, l_group) %>%
  full_join(data_prep %>% 
              expand(nesting(year, stage, l_group))) %>%
  arrange(year, stage)


lam_fit %>% 
  ggplot(aes(year, median, color = stage))+
  geom_line()+
  geom_ribbon(aes(ymin = lower68, ymax = upper68, fill = stage), 
              alpha = 0.2, linetype = 0)+
  scale_color_manual(values=c("dodgerblue2","firebrick2"," black"))+
  scale_fill_manual(values=c("dodgerblue2","firebrick2", "black"))+
  scale_y_continuous(trans = "log", breaks = c(40,80,160,320))


lam_fit %>% 
  ggplot(aes(year, median, color = stage))+
  geom_line()+
  geom_ribbon(aes(ymin = lower68, ymax = upper68, fill = stage),
              alpha = 0.2, linetype = 0)+
  geom_point(data = data_prep %>%
               group_by(year, stage) %>%
               summarize(count = sum(count)),
             aes(y = count))+
  geom_line(data = data_prep %>%
              group_by(year, stage) %>%
              summarize(count = sum(count)),
            aes(y = count),
            size = 0.2, alpha = 0.5)+
  scale_color_manual(values=c("dodgerblue2","firebrick2"," black"))+
  scale_fill_manual(values=c("dodgerblue2","firebrick2", "black"))+
  scale_y_continuous(trans = "log")










nu_fit <- fit_sum %>%
  filter(str_detect(.$rowname, "nu\\[")) %>%
  mutate(n_group = strsplit(rowname, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(mean, sd, lower68, median, upper68, n_group) %>%
  full_join(data_prep %>% 
              expand(nesting(year, stage, n_group))) %>%
  arrange(year, stage)


nu_fit %>% 
  ggplot(aes(year, median, color = stage))+
  geom_line()+
  geom_ribbon(aes(ymin = lower68, ymax = upper68, fill = stage), 
              alpha = 0.2, linetype = 0)+
  scale_color_manual(values=c("dodgerblue2","firebrick2"," black"))+
  scale_fill_manual(values=c("dodgerblue2","firebrick2", "black"))+
  scale_y_continuous(trans = "log", breaks = c(40,80,160,320))


nu_fit %>% 
  ggplot(aes(year, median, color = stage))+
  geom_line()+
  geom_ribbon(aes(ymin = lower68, ymax = upper68, fill = stage), 
              alpha = 0.2, linetype = 0)+
  geom_point(data = data_prep,
             aes(y = count))+
  geom_line(data = data_prep %>%
              group_by(year, stage) %>%
              summarize(count = mean(count)),
            aes(y = count),
            size = 0.2, alpha = 0.5)+
  scale_color_manual(values=c("dodgerblue2","firebrick2"," black"))+
  scale_fill_manual(values=c("dodgerblue2","firebrick2", "black"))

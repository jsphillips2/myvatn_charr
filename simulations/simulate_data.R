#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(nlme)
library(hdrcde)

# import data
model = "char_model"
obs_data = read_csv("data/myvatn_char_counts.csv")  %>%
  mutate(stage = factor(stage, levels=c("first","second","third","adult")))
model_fit = read_csv(paste0("model/output/fit_summary.csv"))
data_list = read_rdump("model/data_list.R")
data_list$sd_obs = data_list$n_obs[2:4,] %>% sd()
fixed_pars = read_csv(paste0("model/output/fixed_pars.csv")) 

# format demographic data
dem_data = model_fit %>% 
  mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult")),
         year = time + 1985) %>%
  filter(name %in% c("phi")) %>%
  select(year,stage,middle) %>%
  spread(stage, middle) %>%
  left_join(model_fit %>% 
              mutate(year = time + 1985) %>%
              filter(name %in% c("rho_scale")) %>%
              select(year,middle) %>%
              rename(rho_scale = middle))

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
#========== Simulation: Fixed Parameters
#==========

# simulation function for fixed parameters
sim_fix_fun = function(data_list_, model_fit_, seed = 1) {
  
  # fixed projection 
  proj = model_fit_ %>%
    filter(name %in% c("phi")) %>%
    mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult"))) %>%
    select(stage, time, middle) %>%
    spread(stage, middle) %>%
    full_join(model_fit_ %>%
                filter(name %in% c("rho_scale")) %>%
                rename(rho_scale = middle) %>%
                select(time, rho_scale)) %>%
    select(-time) %>%
    summarize_all(mean)
  
  # matrix to store data
  mat = matrix(0, nrow = data_list_$nStages, ncol = data_list_$nYears)
  
  # initial population size
  mat[,1] = {model_fit_ %>% 
      filter(name %in% c("N")) %>%
      mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult"))) %>%
      filter(time == min(time,na.rm=T)) %>% 
      mutate(stage = factor(stage, levels=c("first","second","third","adult"))) %>%
      arrange(stage)}$middle
  
  # project demogrpahics
  for (t in 2:ncol(mat)){
    mat[1,t] = proj$rho_scale*mat[4,t-1]
    mat[2,t] = mat[1,t-1]
    mat[3,t] = proj$second*mat[2,t-1]
    mat[4,t] = proj$third*mat[3,t-1] + proj$adult*mat[4,t-1]
  }
  
  # add observation error
  set.seed(seed)
  mat_obs = exp(log(mat) + 
                  rnorm(n = nrow(mat)*ncol(mat), mean = 0, 
                        sd = hdr(fixed_pars$sig_obs_sd, prob = 68)$mode))
  
  # return data
  return(list(mat = mat, mat_obs = mat_obs))
  
}



# test plot: multiple simulations
lapply(1:10, function(x) {
  {sim_fix_fun(data_list, model_fit, x)}$mat_obs %>%
    t() %>% 
    as.data.frame %>%
    rename("first" = "V1", "second" = "V2", "third" = "V3", "adult" = "V4") %>%
    mutate(year = 1:length(first) + 1985) %>%
    gather(stage, count, first, second, third, adult) %>%
    mutate(stage = factor(stage, levels=c("first","second","third","adult")),
           sim = x)
}) %>%
  bind_rows() %>% 
  ggplot(aes(year, count, color = stage))+
  facet_wrap(~stage)+
  geom_line(aes(group = sim), alpha = 0.5, size = 0.3)+
  geom_line(data = {sim_fix_fun(data_list, model_fit, 1)}$mat %>%
              t() %>% 
              as.data.frame() %>%
              rename("first" = "V1", "second" = "V2", "third" = "V3", "adult" = "V4") %>%
              mutate(year = 1:length(first) + 1985) %>%
              gather(stage, count, first, second, third, adult) %>%
              mutate(stage = factor(stage, levels=c("first","second","third","adult"))),
            size = 0.7)+
  scale_color_manual(values=c("gray50","dodgerblue2","firebrick2","black"))+
  scale_y_continuous("Abundance",trans="log",breaks=c(12,48,192))

# simulate and export
n_sim = 100
sim_data_fix = lapply(1:n_sim, function(x) {
  {sim_fix_fun(data_list, model_fit, x)}$mat_obs %>%
    t() %>% as.data.frame() %>%
    rename("first" = "V1", "second" = "V2", "third" = "V3", "adult" = "V4") %>%
    mutate(year = 1:length(first) + 1985) %>%
    gather(stage, count, first, second, third, adult) %>%
    mutate(stage = factor(stage, levels=c("first","second","third","adult")),
           sim = x)
}) %>%
  bind_rows()
# write_csv(sim_data_fix, "simulations/input/fixed/sim_counts.csv")





#==========
#========== Simulation: Observations
#==========

# simulation function for fixed parameters
sim_obs_fun = function(data_list_, model_fit_, seed = 1) {
  
  # fixed projection 
  proj = model_fit_ %>%
    filter(name %in% c("phi")) %>%
    mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult"))) %>%
    select(stage, time, middle) %>%
    spread(stage, middle) %>%
    full_join(model_fit_ %>%
                filter(name %in% c("rho_scale")) %>%
                rename(rho_scale = middle) %>%
                select(time, rho_scale)) 
  
  # matrix to store data
  mat = matrix(0, nrow = data_list_$nStages, ncol = data_list_$nYears)
  
  # initial population size
  mat[,1] = {model_fit_ %>% 
      filter(name %in% c("N")) %>%
      mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult"))) %>%
      filter(time == min(time,na.rm=T)) %>% 
      mutate(stage = factor(stage, levels=c("first","second","third","adult"))) %>%
      arrange(stage)}$middle
  
  # project demogrpahics
  for (t in 2:ncol(mat)){
    mat[1,t] = proj$rho_scale[t-1]*mat[4,t-1]
    mat[2,t] = mat[1,t-1]
    mat[3,t] = proj$second[t-1]*mat[2,t-1]
    mat[4,t] = proj$third[t-1]*mat[3,t-1] + proj$adult[t-1]*mat[4,t-1]
  }
  
  # add observation error
  set.seed(seed)
  mat_obs = exp(log(mat) + 
                  rnorm(n = nrow(mat)*ncol(mat), mean = 0, 
                        sd = hdr(fixed_pars$sig_obs_sd, prob = 68)$mode))
  
  # return data
  return(list(mat = mat, mat_obs = mat_obs))
  
}



# test plot: multiple simulations
lapply(1:10, function(x) {
  {sim_obs_fun(data_list, model_fit, x)}$mat_obs %>%
    t() %>% as.data.frame() %>%
    rename("first" = "V1", "second" = "V2", "third" = "V3", "adult" = "V4") %>%
    mutate(year = 1:length(first) + 1985) %>%
    gather(stage, count, first, second, third, adult) %>%
    mutate(stage = factor(stage, levels=c("first","second","third","adult")),
           sim = x)
}) %>%
  bind_rows() %>% 
  ggplot(aes(year, count, color = stage))+
  facet_wrap(~stage)+
  geom_line(aes(group = sim), alpha = 0.5, size = 0.3)+
  geom_line(data = {sim_obs_fun(data_list, model_fit, 1)}$mat %>%
              t() %>% as.data.frame() %>%
              rename("first" = "V1", "second" = "V2", "third" = "V3", "adult" = "V4") %>%
              mutate(year = 1:length(first) + 1985) %>%
              gather(stage, count, first, second, third, adult) %>%
              mutate(stage = factor(stage, levels=c("first","second","third","adult"))),
            size = 0.7)+
  scale_color_manual(values=c("gray50","dodgerblue2","firebrick2","black"))+
  scale_y_continuous("Abundance",trans="log",breaks=c(12,48,192))

# simulate and export
n_sim = 100
sim_data_obs = lapply(1:n_sim, function(x) {
  {sim_obs_fun(data_list, model_fit, x)}$mat_obs %>%
    t() %>% as.data.frame() %>%
    rename("first" = "V1", "second" = "V2", "third" = "V3", "adult" = "V4") %>%
    mutate(year = 1:length(first) + 1985) %>%
    gather(stage, count, first, second, third, adult) %>%
    mutate(stage = factor(stage, levels=c("first","second","third","adult")),
           sim = x)
}) %>%
  bind_rows()
# write_csv(sim_data_obs, "simulations/input/observation/sim_counts.csv")





#==========
#========== Simulation: Varying Parameters (AR1)
#==========

# function to fit ARMA model with trend
arma_fit_fun = function(data_, par_, struct_, pq_ = NULL) {
  
  # define formula
  form = formula(paste0(par_,"~ time"))
  
  # if ARMA, test different orders and select best
  if (struct_=="ARMA") {
    mods = lapply(1:pq_, function(x){
      m = gls(form, correlation = corARMA(form = ~ time, p = x, q = x),
              data = data_, method = "ML")
      return(m)
    })
    # select model with lowest AIC
    mod_best = mods[which.min(sapply(mods, AIC))][[1]]
  }
  
  # if AR1, fit AR1
  if (struct_=="AR1")  {
    mod_best = gls(form, correlation = corAR1(form = ~ time),
                   data = data_, method = "ML")
  }
  
  # return
  return(mod_best)
}

# function to simulate data from model fit
arma_sim_fun = function(data_, par_) {
  
  # define formula
  form = formula(paste0(par_,"~ time"))
  
  # fit model
  model = gls(form, correlation = corAR1(form = ~ time),
              data = data_, method = "ML")
  
  # extra p's and q's
  phi_theta = coef(model$modelStruct$corStruct, unconstrained = F)
  
  
  # construct ARMA structure
  arma_struct = list(order = c(1, 0, 0), ar = phi_theta)
  
  # simualte ARMA without trend
  sims = arima.sim(n = length(model$fitted), model = arma_struct,
                   rand.gen =  function(n, ...) rnorm(n, mean = 0, sd = model$sigma))
  
  # add trend
  sims_trend = c(sims + predict(model))
  
  # return simulated data
  return(sims_trend)
  
}



# define invserse logit function
inv_logit = function(x){exp(x)/(1 + exp(x))}



# simulation function for fixed parameters
sim_var_fun = function(data_list_, model_fit_, seed = 1) {
  
  # set seed
  set.seed(seed)
  
  # fixed projection 
  proj = model_fit %>%
    filter(name %in% c("logit_phi")) %>%
    mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult"))) %>%
    select(stage, time, middle) %>%
    spread(stage, middle) %>%
    full_join(model_fit %>%
                filter(name %in% c("log_rho")) %>%
                rename(log_rho = middle) %>%
                select(time, log_rho)) %>%
    arrange(time)
  
  # calculate rho scale
  scale = {model_fit %>% 
      filter(name %in% c("rho","rho_scale")) %>% 
      select(name,time, middle) %>% 
      spread(name, middle) %>% 
      mutate(scale = rho_scale/rho)}$scale[1]
  
  # simulate demographic parameters
  proj_sim = tibble(
    second = arma_sim_fun(proj, "second") %>% inv_logit(),
    third = arma_sim_fun(proj, "third")  %>% inv_logit(),
    adult = arma_sim_fun(proj, "adult")  %>% inv_logit(),
    rho_scale = scale*{arma_sim_fun(proj, "log_rho") %>% exp()}
  )
  
  # matrix to store data
  mat = matrix(0, nrow = data_list_$nStages, ncol = data_list_$nYears)
  
  # initial population size
  mat[,1] = {model_fit_ %>% 
      filter(name %in% c("N")) %>%
      mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult"))) %>%
      filter(time == min(time,na.rm=T)) %>% 
      mutate(stage = factor(stage, levels=c("first","second","third","adult"))) %>%
      arrange(stage)}$middle
  
  # project demogrpahics
  for (t in 2:ncol(mat)){
    mat[1,t] = proj_sim$rho_scale[t-1]*mat[4,t-1]
    mat[2,t] = mat[1,t-1]
    mat[3,t] = proj_sim$second[t-1]*mat[2,t-1]
    mat[4,t] = proj_sim$third[t-1]*mat[3,t-1] + proj_sim$adult[t-1]*mat[4,t-1]
  }
  
  # add observation error
  mat_obs = exp(log(mat) + 
                  rnorm(n = nrow(mat)*ncol(mat), mean = 0, 
                        sd = hdr(fixed_pars$sig_obs_sd, prob = 68)$mode))
  
  # return data
  return(mat_obs)
  
}



# test plot: multiple simulations
lapply(1:20, function(x) {
  {sim_var_fun(data_list, model_fit, x)} %>%
    t() %>% as.data.frame() %>%
    rename("first" = "V1", "second" = "V2", "third" = "V3", "adult" = "V4") %>%
    mutate(year = 1:length(first) + 1985) %>%
    gather(stage, count, first, second, third, adult) %>%
    mutate(stage = factor(stage, levels=c("first","second","third","adult")),
           sim = x)
}) %>%
  bind_rows() %>% 
  ggplot(aes(year, count, color = stage))+
  facet_wrap(~stage)+
  geom_line(aes(group = sim), alpha = 0.5, size = 0.3)+
  scale_color_manual(values=c("gray50","dodgerblue2","firebrick2","black"))+
  scale_y_continuous("Abundance",trans="log",breaks=c(12,48,192))

# simulate and export
n_sim = 100
sim_data_var = lapply(1:n_sim, function(x) {
  {sim_var_fun(data_list, model_fit, x)} %>%
    t() %>% as.data.frame() %>%
    rename("first" = "V1", "second" = "V2", "third" = "V3", "adult" = "V4") %>%
    mutate(year = 1:length(first) + 1985) %>%
    gather(stage, count, first, second, third, adult) %>%
    mutate(stage = factor(stage, levels=c("first","second","third","adult")),
           sim = x)
}) %>%
  bind_rows()
# write_csv(sim_data_var, "simulations/input/variable/sim_counts.csv")






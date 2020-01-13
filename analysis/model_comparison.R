#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(loo)

# read data
data <- read_rdump("model/data_list.R")





#==========
#========== LOO
#==========

# model names
models <- c("model_full","model_rho","model_phi","model_fixed") 

# loo
loos <- lapply(models, function(x){

  fit <- readRDS(paste0("model/output/",x,"/fit.rds"))

  log_lik <- extract_log_lik(fit, merge_chains = FALSE)
  r_eff <- relative_eff(exp(log_lik))

  loo <- loo(log_lik, r_eff = r_eff, cores = 10)
}) %>%
  set_names(models)

# looic
looic <- compare(loos[[1]],
                 loos[[2]],
                 loos[[3]],
                 loos[[4]]) %>%
  as.data.frame() %>%
  rownames_to_column(var = "model")  %>%
  as_tibble() %>%
  # select(model, elpd_loo, elpd_diff, p_loo, looic) %>%
  mutate(model = models[str_split(model, "\\[\\[|\\]\\]") %>% 
                          map_int(~as.integer(.x[2]))]) %>%
  mutate(delt_looic = looic - looic[1])
  




#==========
#========== Log likelihoods
#==========


# model names
models <- c("model_full","model_rho","model_phi","model_fixed") 

# log likelihood
log_liks <- lapply(models, function(x){
  
  fit <- readRDS(paste0("model/output/",x,"/fit.rds"))
  
  chains <- fit@stan_args %>% length()
  iter <- fit@stan_args[[1]]$iter
  fit_summary <- summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
  {as_tibble(.) %>%
      mutate(var = rownames(summary(fit)$summary))}
  
  log_liks <- tibble(model = x,
                     mean_log_lik = {fit_summary %>%
                         filter(var=="log_lik_sum")}$`50%`
                     )
  return(log_liks)
  
}) %>%
  bind_rows()





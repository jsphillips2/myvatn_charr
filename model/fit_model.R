#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# read data
data <- read_rdump("model/data_list.R")

# introduce systematic bias in age 2 fish
bias <- F
if(bias == T){data$N_obs <- exp(data$n_obs)
data$N_obs[2,] <- 2*data$N_obs[2,]
data$n_obs <- log(data$N_obs)}

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
#========== Fit model
#==========

# select model
models <- c("model_full","model_rho","model_phi","model_fixed") 
model <- models[1]
model_path <- paste0("model/stan/",model,".stan")
export_path <- if(bias == T){paste0("model/output/",model,"_bias")
} else{paste0("model/output/",model)}

# MCMC specificaions
chains <- 4
iter <- 2000 
adapt_delta <- 0.975
max_treedepth <- 15

# fit model
fit <- stan(file = model_path, data = data, seed=1, chains = chains, iter = iter, 
            control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))

# summary of fit
fit_summary <- summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>% 
{as_tibble(.) %>%
    mutate(var = rownames(summary(fit)$summary))}

# check Rhat & n_eff
fit_summary %>% filter(Rhat > 1.01) %>% select(Rhat, n_eff, var) %>% arrange(-Rhat)
fit_summary %>% filter(n_eff < 0.5*(chains*iter/2)) %>% select(Rhat, n_eff, var) %>% arrange(n_eff) %>%
  mutate(eff_frac = n_eff/(chains*iter/2))

# save model full output
# saveRDS(fit, paste0(export_path,"/fit.rds"))

# import model fit 
# fit <- readRDS(paste0("model/output/",model,"/fit.rds"))
# chains <- fit@stan_args %>% length()
# iter <- fit@stan_args[[1]]$iter
# fit_summary <- summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
# {as_tibble(.) %>%
#     mutate(var = rownames(summary(fit)$summary))}







#==========
#========== Examine Chains
#==========

# function for selecting fixed parameters
fixed_par_fn <- function(x){
  if (x == "model_full"){return(c("log_rho[1]","logit_phi[1,1]","logit_phi[2,1]","logit_phi[3,1]",
                                  "sig_rho","sig_phi","sig_obs","sig_obs_sd",
                                  "lp__"))}
  if (x == "model_rho"){return(c("log_rho[1]","logit_phi[1]","logit_phi[2]","logit_phi[3]",
                                 "sig_rho","sig_obs","sig_obs_sd",
                                 "lp__"))}
  if (x == "model_phi"){return(c("log_rho","logit_phi[1,1]","logit_phi[2,1]","logit_phi[3,1]",
                                 "sig_phi","sig_obs","sig_obs_sd",
                                 "lp__"))}
  if (x == "model_fixed"){return(c("log_rho","logit_phi[1]","logit_phi[2]","logit_phi[3]",
                                   "sig_obs","sig_obs_sd",
                                   "lp__"))}
}

# extract fixed parameters
fixed_pars <- rstan::extract(fit, pars=fixed_par_fn(model)) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(fixed_par_fn(model)) %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))

# examine chains for parameters
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  ggplot(aes(step, value, color=factor(chain)))+
  facet_wrap(~par, scales="free_y")+
  geom_line(alpha=0.5)+
  theme_bw()

# pairs plot for parameters
# GGally::ggpairs(fixed_pars %>% select(-chain, -step, -lp__, -sig_obs_sd))

# pairs plot for parameters
GGally::ggpairs(fixed_pars %>% select(sig_rho, sig_phi, sig_obs))


# posterior densities
fixed_pars %>%
  gather(par, value, -chain, -step) %>%
  filter(par != "lp__") %>%
  ggplot(aes(value))+
  facet_wrap(~par, scales="free")+
  stat_density(alpha=0.5, geom = "line")+
  theme_bw()






#==========
#========== Posterior Predictive Check (residual autocorrelation)
#==========

# select autocorrelations
post_pars <- c(crossing(stage = 1:3) %>%
{paste0("acor_obs[",.$stage,"]")},
crossing(stage = 1:3) %>%
{paste0("acor_sim[",.$stage,"]")}) 

# extract autocorrelations
post_pred = rstan::extract(fit, pars=post_pars) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(post_pars) %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains))

# plot
post_pred %>%
  gather(var, val, -chain, -step) %>%
  mutate(type = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         stage = strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2]))) %>%
  select(-var) %>%
  spread(type, val) %>%
  ggplot(aes(acor_obs, acor_sim))+
  facet_wrap(~stage, nrow = 3)+
  geom_point(alpha = 0.5)+
  # geom_abline(intercept = 0, slope = 1)+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  # scale_y_continuous("Residual Autocorrelation (Simulated)",limits=c(-0.7,0.7))+
  # scale_x_continuous("Residual Autocorrelation (Observed)",limits=c(-0.7,0.7))+
  coord_equal()





#==========
#========== Clean output and export
#==========

# select vars to extract for demographic analysis
time_pars_fn <- function(x){
  if(x == "model_full"){return(c(paste0("rho_scale[",1:(data$nYears-1),"]"),
                                 crossing(stage = 1:3, time = 1:(data$nYears-1)) %>%
                                 {paste0("phi[",.$stage,",",.$time,"]")}))}
  if(x == "model_rho"){return(c(paste0("rho_scale[",1:(data$nYears-1),"]"),
                              crossing(stage = 1:3) %>%
                                {paste0("phi[",.$stage,"]")}))}
  if(x == "model_phi"){return(c("rho_scale",
                                crossing(stage = 1:3, time = 1:(data$nYears-1)) %>%
                                 {paste0("phi[",.$stage,",",.$time,"]")}))}
  if(x == "model_fixed"){return(c("rho_scale",
                                  crossing(stage = 1:3) %>%
                                  {paste0("phi[",.$stage,"]")}))}
}

# extract full chains for demographic analysis
time_pars <- rstan::extract(fit, pars=time_pars_fn(model)) %>%
  lapply(as_tibble) %>%
  bind_cols() %>%
  set_names(time_pars_fn(model)) %>%
  mutate(chain = rep(1:chains, each = iter/2), step = rep(c(1:(iter/2)), chains)) %>%
  gather(var, value, -chain, -step) %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         stage = ifelse(name %in% c("phi"),
                        strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])), NA),
         time = ifelse(name %in% c("phi"),
                       strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
                       ifelse(name %in% c("rho_scale"),
                              strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])), NA))) %>%
  select(name, chain, step, stage, time, value) %>%
  arrange(name, chain, step, stage, time)

# summary
fit_clean <- fit_summary %>%
  rename(lower16 = `16%`, middle = `50%`, upper84 = `84%`)  %>%
  mutate(name = strsplit(var, "\\[|\\]|,") %>% map_chr(~.x[1]),
         stage = ifelse(name %in% c("n_init","n_pred","N"),
                        strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])), 
                        ifelse(name %in% c("logit_phi_init","logit_phi","phi","sig_phi"),
                               strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])+1L), NA)),
         time = ifelse(name %in% c("logit_phi_init","logit_phi","phi","n_init","n_pred","N"),
                       strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])), 
                       ifelse(name %in% c("log_rho","rho","N_tot","lambda","rho_scale"),
                              strsplit(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])), NA))) %>%
  select(name, stage, time, middle, lower16, upper84) 

# export
export_path <- if(bias == T){paste0("model/output/",model,"_bias")
  } else{paste0("model/output/",model)}
# write_csv(fixed_pars, paste0(export_path, "/fixed_pars.csv"))
# write_csv(fit_clean, paste0(export_path,"/fit_summary.csv"))
# write_csv(time_pars, paste0(export_path,"/time_pars.csv"))
# write_csv(post_pred, paste0(export_path,"/post_pred_pars.csv"))










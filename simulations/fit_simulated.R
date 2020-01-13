#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(hdrcde)

# stan settings
rstan_options(auto_write = TRUE)

# simulation type
type = "fixed"

# import data
data_full = read_csv(paste0("simulations/input/",type,"/sim_counts.csv"))

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
#========== Functions for fitting simulated data
#==========

# function to prepare data for fit
prep_fun = function(data_, x_) {
  
  #  select data for simulation x
  d = data_ %>% 
    filter(sim == x_)
  
  # spread data
  d_wide = d %>%
    # log transform
    mutate(count = log(count)) %>%
    # convert to wide format
    spread(stage, count) 
  
  # define list elements
  n_obs = d_wide[,c("first","second","third","adult")] %>% t()
  nYears = ncol(n_obs)
  nStages = nrow(n_obs)
  
  # data list
  d_list = list(nYears = nYears, nStages = nStages, n_obs = n_obs)
  
  # return list
  return(d_list)
}



# function for fitting simulated data
fit_fun = function(data_, model_ = "char_model", chains_ = 1, iter_ = 1000, cores = 1) {
  
  # define model path
  model_path = paste0("model/",model_,".stan")
  
  # fit model
  fit = stan(file = model_path, data = data_, chains = chains_, iter = iter_, 
             control = list(adapt_delta = 0.9))
  
  # return fit
  return(fit)
  
}

# function for checking divergences
check_div_custom <- function(fit_) {
  
  # extract divergences
  sampler_params <- get_sampler_params(fit_, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  # return
  return(n/N)
}

# function for processing output
proc_fun = function(fit_) {
  
  # extract relevant parameters
  fixed_par_v = c("sig_rho","sig_phi","sig_obs","sig_obs_sd","lp__")
  fixed_pars = rstan::extract(fit_, pars=fixed_par_v) %>%
    lapply(as_tibble) %>%
    bind_cols() %>%
    set_names(fixed_par_v) 
  
  # posterior summary
  fit_clean = fixed_pars %>%
    gather(name, val, sig_rho, sig_phi, sig_obs, sig_obs_sd) %>%
    group_by(name) %>%
    summarize(lower16 = hdr(val, prob = 68)$hdr[1],
              middle = hdr(val, prob = 68)$mode,
              upper84 = hdr(val, prob = 68)$hdr[2])
  
  # rhat
  fit_summary = summary(fit_)$summary 
  fit_summary = {as_tibble(fit_summary) %>%
      mutate(var = rownames(fit_summary))} %>%
    filter(var %in% c("sig_rho","sig_phi","sig_obs","sig_obs_sd"))
  Rhat_max = max(fit_summary$Rhat)
  
  # divergence
  div = check_div_custom(fit_)
  
  
  # return clean fit
  return(fit_clean %>% mutate(Rhat_max = Rhat_max, div = div))
}





#==========
#========== Fit to simulated data
#==========

# set number of cores 
options(mc.cores = parallel::detectCores()-2)

# fit simulated data
set.seed(1)
fit_sim1 = parallel::mclapply(1:50, function(x) {
  
  # fit and summarize data for simulation x
  d = proc_fun(fit_fun(prep_fun(data_full, x))) %>%
    mutate(sim = x)
  
  # return fit
  return(d)
  
}) %>%
  bind_rows()

fit_sim2 = parallel::mclapply(51:100, function(x) {
  
  # fit and summarize data for simulation x
  d = proc_fun(fit_fun(prep_fun(data_full, x))) %>%
    mutate(sim = x)
  
  # return fit
  return(d)
  
}) %>%
  bind_rows()

fit_sim = bind_rows(fit_sim1, fit_sim2)

# export
# write_csv(fit_sim, paste0("simulations/output/", type, "/fit_sim.csv"))

# calculate density function weights by se of estimate
fit_density = fit_sim %>%
  split(.$name) %>%
  lapply(function(x){
    
    # calcualte density
    d = density(x = x$middle, from = 0)
    
    # return
    return(tibble(name = x$name %>% unique(), est = d$x, dens = d$y))
  }) %>%
  bind_rows()

# extract original paramters estimates
orig_est = read_csv(paste0("model/output/fixed_pars.csv")) %>%
  gather(name, val, sig_rho, sig_phi, sig_obs, sig_obs_sd) %>%
  group_by(name) %>%
  summarize(lower16 = hdr(val, prob = 68)$hdr[1],
            middle = hdr(val, prob = 68)$mode,
            upper84 = hdr(val, prob = 68)$hdr[2]) %>%
  mutate(ymin = 0) %>%
  left_join(fit_density %>%
              group_by(name) %>%
              summarize(ymax = max(dens)) %>% 
              ungroup())

# plot
sd_exc = "sig_obs"
fit_density %>%
  filter(name != sd_exc) %>%
  ggplot()+
  facet_wrap(~name, scales="free", nrow = 3)+
  geom_line(aes(est, dens))+
  geom_segment(data  = orig_est %>% filter(name != sd_exc), 
               aes(x = middle, xend = middle, y = ymin, yend = ymax), 
               size = 0.6, linetype = 2)+
  geom_rect(data = orig_est %>% filter(name != sd_exc), 
            aes(xmin = lower16, xmax = upper84, ymin = ymin, ymax = ymax),
            alpha = 0.2)+
  scale_y_continuous("Probability Density")+
  scale_x_continuous("Estimate")

# p-value
{fit_density %>%
  filter(name =="sig_obs_sd") %>%
  mutate(pp = ifelse(est > {orig_est %>% filter(name=="sig_obs_sd")}$middle, 1, 0))}$pp %>% mean()


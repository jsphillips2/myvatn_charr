#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(rstan)
library(Matrix)

# stan settings
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

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
#========== Prepare Data
#==========

data_prep = data %>%
  #Fill gaps in "Est.age"
  mutate(est_age = ifelse(is.na(age)==F, age, est_age)
  ) %>%
  #Remove June data & samples without ages
  filter(month != 6, is.na(est_age)==F) %>%
  #Define age_group
  mutate(age_group = factor(ifelse(est_age == 1, "first", 
                               ifelse(est_age==2, "second", 
                                      ifelse(est_age==3, "third", "adult"))),
                        levels = c("first","second","third","adult"))
  )

# trim data
data_trim <- data_prep %>%
  filter(!(is.na(mass)),
         mass > 0,
         !(mass > 1000 & age_group == "second"))

# mean mass through time, by age
data_trim %>%
  group_by(year, age_group) %>%
  summarize(lo = quantile(mass, probs = 0.16),
            hi = quantile(mass, probs = 0.84),
            mass = median(mass)) %>%
  ggplot(aes(year, mass, color = age_group))+
  geom_line()+
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = age_group), linetype = 0, alpha = 0.2) +
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2"))+
  scale_fill_manual(values = c("gray50","goldenrod", "firebrick","skyblue2"))+
  scale_y_continuous(Mass~(g),
                     trans = "log", breaks = c(50,200,800,2400))

# mass distribution, by age
data_trim %>%
  filter(mass > 6) %>%
  ggplot(aes(mass, fill = age_group))+
  geom_vline(xintercept = 35)+
  geom_histogram(position = "identity",alpha = 0.5, binwidth = 0.1)+
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2"))+
  scale_fill_manual(values = c("gray50","goldenrod", "firebrick","skyblue2"))+
  scale_x_continuous(trans = "log",  breaks = c(50,200,800,2400))

# mass at age, through time
data_trim %>%
  filter(mass > 6) %>%
  mutate(group = ifelse(year > 2000, "b", "a")) %>%
  ggplot(aes(est_age, mass, color = year))+
  geom_hline(yintercept = 400)+
  geom_jitter(width = 0.2, alpha = 0.5, size= 1)+
  scale_y_continuous(Mass~(g),
                     trans = "log", breaks = c(50,200,800,2400))+
  scale_color_gradient2(low = "firebrick", 
                        high = "dodgerblue", 
                        mid = "gray50", 
                        midpoint = 2000)

# mass quantiles (assume minimum catch mass of 35)
nsize = 3
bins <- seq(0, 1, length.out = nsize + 1)
mqq <- {data_trim %>% filter(mass > 35)}$mass %>% quantile(probs = bins)
mq <- c(min(data_trim$mass), mqq)  
mp <- (mq %>% diff())/2 + mq[1:(nsize + 1)] %>% round(1)

# plot
data_trim %>%
  filter(mass > 6) %>%
  ggplot(aes(mass, fill = age_group))+
  geom_histogram(position = "identity",alpha = 0.5, binwidth = 0.1)+
  geom_vline(data = tibble(x = mq[2:(nsize+1)]), aes(xintercept = x))+
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2"))+
  scale_fill_manual(values = c("gray50","goldenrod", "firebrick","skyblue2"))+
  scale_x_continuous(trans = "log",  breaks = as.numeric(mq[2:(nsize+1)]))

# count data by age x mass
counts <- data_trim %>%
  mutate(mass_class = cut(data_trim$mass, as.numeric(mq)) %>% as.numeric(),
         age_class = cut(est_age, c(0,1,2,3,max(est_age))) %>% as.numeric()) %>%
  select(year, age_class, mass_class, mass) %>%
  na.omit() %>%
  group_by(year, age_class, mass_class) %>%
  summarize(count = length(mass)) %>%
  ungroup()

# expand count data
counts_expand <- counts %>%
  full_join(tidyr::crossing(year = seq(min(counts$year), max(counts$year), by = 1),
                          mass_class = seq(min(counts$mass_class), max(counts$mass_class), by = 1), 
                          age_class = seq(min(counts$age_class), max(counts$age_class), by = 1)))  %>%
  mutate(count = ifelse(is.na(count)==T, 0L, count)) %>%
  arrange(year, age_class, mass_class)

# plot
counts_expand %>%
  ggplot(aes(year, count, color = factor(age_class), linetype = factor(mass_class)))+
  geom_line()+
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2","black"))





#==========
#========== Package Data
#==========


### General set up
masses <- 4 # number of mass classes
ages <- 4 # number  of age classes
times <- counts_expand$year %>% unique() %>% length() # number of time steps
q <- mp/mean(mp) # mass covariate (scale by mean)
fa <- 4 # minimum reproductive age



### Observed data
y <- matrix(counts_expand$count, 
            nrow = ages * masses, 
            ncol = times)



### Observation status
yo <- {counts_expand %>%
    select(age_class, mass_class) %>%
    unique() %>%
    mutate(obs = ifelse(mass_class == 1, 0, 1))}$obs


# Priors
p <- matrix(c(1.5,     1.5/0.5,      # b0s
              1.5,     1.5/0.5,      # b1s
              1.5,     1.5/0.5,      # us
              1.5,     1.5/0.5,      # b0
              1.5,     1.5/0.5,      # b1
              1.5,     1.5/0.5,      # u
              1.5,     1.5/0.5,      # f
              1.5,     1.5/10     # ys
), nrow = 8, ncol = 2, byrow = T)

# weigting for penalization
ws <- 4


#### Matrices
# function for block-diagonal matrix from submatrix
block_fn <- function(v_, r_, c_) {
  as.matrix(Reduce("bdiag", 
                   lapply(1:c_, function(x){
                     matrix(v_, nrow = r_, ncol = r_, byrow = T)
                   })))
}

D <- block_fn(v_ = c(0, 0, 0, 0,
                     1, 0, 0, 0,
                     0, 1, 0, 0,
                     0, 0, 1, 1),
                    c_ = masses, r_ = ages)

H <- block_fn(v_ = c(1, 1, 1, 1,
                     0, 0, 0, 0,
                     0, 0, 0, 0,
                     0, 0, 0, 0),
                    c_ = masses, r_ = ages)


# create vec-perumtation matrix
k_fn <- function(i_, j_, masses_, ages_) {
  eij_ = matrix(0, nrow = masses_, ncol = ages_)
  eij_[i_, j_] = 1
  return(kronecker(eij_, t(eij_)))
}

ijs <- expand.grid(1:masses, 1:ages)

K = Reduce('+', 
                 lapply(1:nrow(ijs), function(x){
                   k_fn(ijs[x,1], ijs[x,2], masses, ages)
                 }))

# package data
data_list <- list(times = times,
                  masses = masses,
                  ages = ages,
                  q = q,
                  fa = fa,
                  K = K,
                  D = D,
                  H = H,
                  y = y,
                  yo = yo,
                  p = p,
                  ws = ws)





#==========
#========== Fit model
#==========

# models
models <- c("size_age_model")
model <- models[1]

# model
model_path <- paste0("multistate_model/",model,".stan")

# init_fn <- function() {
#   list(b0s = runif(1, 0.5 * p[1, 1]/p[1, 2], 1.5 * p[1, 1]/p[1, 2]),
#        b1s = runif(1, 0.5 * p[2, 1]/p[2, 2], 1.5 * p[2, 1]/p[2, 2]),
#        us = runif(1, 0.5 * p[3, 1]/p[3, 2], 1.5 * p[3, 1]/p[3, 2]),
#        b0 = runif(1, 0.5 * p[4, 1]/p[4, 2], 1.5 * p[4, 1]/p[4, 2]),
#        b1 = runif(1, 0.5 * p[5, 1]/p[5, 2], 1.5 * p[5, 1]/p[5, 2]),
#        u = matrix(runif((masses - 1) * (times - 1),
#                          0.5 * p[6, 1]/p[6, 2], 1.5 * p[6, 1]/p[6, 2]),
#                    nrow = masses - 1, ncol = times - 1),
#        f = runif(1, 0.5 * p[7, 1]/p[7, 2], 1.5 * p[7, 1]/p[7, 2]),
#        ys = runif(1, 0.5 * p[8, 1]/p[8, 2], 1.5 * p[8, 1]/p[8, 2]),
#        x0  = runif(ages * masses, 0.5 * mean(y), 1.5 * mean(y)))
# }


# MCMC specifications (for testing)
# chains <- 1
# iter <- 50
# adapt_delta <- 0.8
# max_treedepth <- 10

# MCMC specifications
# chains <- 4
# iter <- 500
# adapt_delta <- 0.9
# max_treedepth <- 10

# fit model
# fit <- stan(file = model_path, data = data_list, seed=2e3,
#             chains = chains, iter = iter,
#             control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth))

# saveRDS(fit, paste0("multistate_model/",model,"_fit.RDS"))

# fit <- readRDS("multistate_model/size_age_model_fixed_fit.RDS")
fit_summary <- rstan::summary(fit, probs=c(0.16, 0.5, 0.84))$summary %>%
{as_tibble(.) %>%
    mutate(var = rownames(rstan::summary(fit)$summary))}



x <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "x"),
         !(str_detect(fit_summary$var, "x0"))) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3]))) %>%
  left_join(counts_expand %>%
              select(age_class, mass_class) %>%
              unique() %>%
              mutate(id = row_number())) %>%
  left_join(tibble(year = counts_expand$year %>% unique(),
                   time = 1:times))


x %>%
  ggplot(aes(year, mi, color = factor(age_class), linetype = factor(mass_class)))+
  facet_wrap(~age_class, scales = "free_y")+
  geom_line()+
  geom_point(data = counts_expand, aes(y = count))+
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2","black"))



x %>%
  ggplot(aes(year, mi, color = factor(age_class), linetype = factor(mass_class)))+
  facet_wrap(~age_class, scales = "free_y")+
  geom_line()+
  geom_ribbon(aes(ymin = lo, ymax = hi, 
                  fill = factor(age_class), 
                  group = factor(mass_class)),
              linetype = 0, alpha = 0.2)+
  geom_point(data = counts_expand, aes(y = count))+
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2","black"))+
  scale_fill_manual(values = c("gray50","goldenrod", "firebrick","skyblue2","black"))



AA <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "AA")) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         row = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         col = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[4]))) %>%
  left_join(counts_expand %>%
              select(age_class, mass_class) %>%
              unique() %>%
              mutate(col = row_number()) %>%
              rename(donor_age = age_class,
                     donor_mass = mass_class)) %>%
  left_join(counts_expand %>%
              select(age_class, mass_class) %>%
              unique() %>%
              mutate(row = row_number()) %>%
              rename(recipient_age = age_class,
                     recipient_mass = mass_class)) %>%
  left_join(tibble(year = counts_expand$year %>% unique(),
                   time = 1:times))



AA %>%
  filter(donor_age == 4) %>%
  ggplot(aes(year, mi, 
             color = factor(recipient_age), 
             linetype = factor(recipient_mass),
             group = interaction(recipient_mass, donor_mass)))+
  facet_wrap(~recipient_age)+
  geom_line()+
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2","black"))



x %>%
  group_by(year, age_class) %>%
  summarize(count = sum(mi)) %>%
  ungroup() %>%
  mutate(cohort = year - age_class,
         cohort = cohort - min(cohort) + 1) %>%
  ggplot(aes(age_class, count, group = cohort, color = year))+
  geom_line()+ 
  geom_point()



x %>%
  group_by(year, age_class) %>%
  summarize(count = sum(mi)) %>%
  ungroup() %>%
  ggplot(aes(year, count, color = factor(age_class)))+
  geom_point(data = counts_expand %>%
               group_by(year, age_class) %>%
               summarize(count = sum(count)))+
  geom_line(data = counts_expand %>%
               group_by(year, age_class) %>%
               summarize(count = sum(count)),
            size = 0.2)+
  geom_line(size = 1)+ 
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2","black"))


fit_summary %>% 
  filter(var == "f")

u <- fit_summary %>%
  select(var, `16%`, `50%`, `84%`) %>%
  rename(lo = `16%`, mi = `50%`, hi = `84%`) %>%
  filter(str_detect(fit_summary$var, "u"),
         !(str_detect(fit_summary$var, "us"))) %>%
  mutate(name = str_split(var, "\\[|\\]|,") %>% map_chr(~as.character(.x[1])),
         id = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[2])),
         time = str_split(var, "\\[|\\]|,") %>% map_int(~as.integer(.x[3]))) %>%
  left_join(counts_expand %>%
              select(mass_class) %>%
              unique() %>%
              mutate(id = row_number())) %>%
  left_join(tibble(year = counts_expand$year %>% unique(),
                   time = 1:times)) %>%
  na.omit()


u %>%
  ggplot(aes(year, mi, color = factor(mass_class)))+
  geom_line()+
  geom_ribbon(aes(ymin = lo, ymax = hi,
                  fill = factor(mass_class)),
              linetype = 0, alpha = 0.2)+
  scale_color_manual(values = c("gray50","goldenrod", "firebrick","skyblue2","black"))+
  scale_fill_manual(values = c("gray50","goldenrod", "firebrick","skyblue2","black"))


bayesplot::mcmc_trace(x = fit, pars = "b1")

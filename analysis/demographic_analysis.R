#==========
#========== Preliminaries
#==========

# load packages
library(tidyverse)
library(hdrcde)
library(matrixcalc)

# set number of cores
options(mc.cores = parallel::detectCores()-2)

# specify model
models <- c("model_full","model_rho","model_phi","model_full_bias") 
model <- models[1]

# import data
time_pars = read_csv(paste0("model/output/",model,"/time_pars.csv"),
                     col_type = "cdddid") %>%
  mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult")),
         year = time + 1985) 
model_fit = read_csv(paste0("model/output/",model,"/fit_summary.csv")) %>%
  mutate(stage = factor(stage, levels = 1:4, labels=c("first","second","third","adult")),
         year = time + 1985) 

# extract demogrpahic pars
set.seed(1)
dem_pars = time_pars  %>%
  filter(name %in% c("phi")) %>%
  spread(stage, value) %>%
  select(-name, -time) %>%
  {if(is.na(.$year[1])==T){select(.,-year)} 
    else return(.)} %>%
  left_join(time_pars %>% 
              filter(name == "rho_scale") %>%
              rename(rec = value) %>%
              select(chain, step, year, rec) %>%
              {if(is.na(.$year[1])==T){select(.,-year)} 
                else return(.)}) %>%
  rename(s32 = first,
         s43 = second,
         s44 = third) %>%
  mutate(s21 = 1,
         chain_step = as.numeric(paste0(chain, step))) %>%
  filter(chain_step %in% sample(chain_step, 1000))





#==========
#========== Define functions: sensitivies, elasticities, and lambda
#==========

# function to calculate elasticities and sensitivies
ela_sens_fn = function(theta){
  
  # parameter vector
  pv = theta %>% 
    select(rec,s21,s32,s43,s44) %>%
    as.numeric
  
  # parameter matrix
  pmat = matrix(rep(pv[2:5],each=4), ncol=4, nrow=4)
  
  # life cycle graph
  cmat = matrix(c(0,1,0,0,0,0,1,0,0,0,0,1,1,0,0,1), ncol=4, nrow=4)
  
  # projection matrix
  mat = hadamard.prod(pmat, cmat)
  mat[1,4] = pv[1]
   
  # extract eigenvalues and eigenvectors
  val = eigen(mat)$values
  W = eigen(mat)$vectors
  imax = which.max(Re(val))
   
  # calculate left eigenvectors
  V = Conj(solve(W))
   
  # select eigenvectors corresponding to max eigenvalue
  w=W[,imax]
  v=Re(V[imax,])
   
  # calculate sensitivies
  smat = v %*% t(w) 
  
  # replace nonsense elements with 0
  ses = hadamard.prod(Re(smat), cmat)
  
  # calculate lambda
  lam = val[imax]
  
  # calculate elasticities
  els = (hadamard.prod(Re(smat), mat)/Re(lam))
  
  # print results
  results = 
    list(sensitivity = ses, elasticity = els, lambda = lam, 
         year = theta$year, chain = theta$chain, step = theta$step)
  
  return(results)
  
}




# function to clean output
ela_out_fn = function(x){
  data.frame(
    year = x$year %>% as.numeric(),
    s_rec = x$sensitivity[1,4],
    s_s21 = x$sensitivity[2,1],
    s_s32 = x$sensitivity[3,2],
    s_s43 = x$sensitivity[4,3],
    s_s44 = x$sensitivity[4,4],
    e_rec = x$elasticity[1,4],
    e_s21 = x$elasticity[2,1],
    e_s32 = x$elasticity[3,2],
    e_s43 = x$elasticity[4,3],
    e_s44 = x$elasticity[4,4],
    lambda = Re(x$lambda)
  )
}





#==========
#========== Calculate sensitivies, elasticities, and lambda
#==========

# calculate sensitivities adn elasticities
sensitivities = dem_pars %>%
  
  # add index for row
  mutate(row = as.numeric(paste0(chain, step))) %>%
  
  # iterate by year
  split(.$year) %>%
  parallel::mclapply(function(x){
    
    # iterate by row
    x %>% split(.$row) %>%
      lapply(function(xx){
        
        # calculate 
        ela_sens_fn(xx) %>%
          ela_out_fn
      }) %>%
      bind_rows() %>%
      gather(name, val, -year) %>%
      group_by(year, name) %>%
      summarize(lower16 = quantile(val, 0.18),
                middle = median(val),
                upper84 = quantile(val, 0.84)) %>%
      ungroup()
    
  }) %>%
  bind_rows() %>%
  mutate(
    type = ifelse(
      name == "e_rec"|name == "e_s21"|name == "e_s32"|
        name == "e_s43"|name == "e_s44", "elasticity",
      ifelse(
        name == "s_rec"|name == "s_s21"|name == "s_s32"|
          name == "s_s43"|name == "s_s44", "sensitivity","lambda"
      )
    )
  )

# export
write_csv(sensitivities, paste0("analysis/demographic_output/",model,"/sensitivities.csv"))



#==========
#========== Asymptotic growth rate
#==========

# calculate reference demographic matrix based on mean values
dem_ref = dem_pars %>% 
  group_by(chain, step) %>%
  summarize(
    rec = mean(rec), s21 = mean(s21), s32 = mean(s32), 
    s43 = mean(s43), s44 = mean(s44)
  ) %>%
  ungroup() %>%
  summarize(
    rec = mean(rec),
    s21 = mean(s21),
    s32 = mean(s32),
    s43 = mean(s43),
    s44 = mean(s44)
  ) %>%
  mutate(chain = NA, step = NA, year = NA)

# add asymptotic growth rate for the mean matrix
dem_ref$lam_reference = Re(ela_sens_fn(dem_ref)$lambda)

# export
write_csv(dem_ref, paste0("analysis/demographic_output/",model,"/dem_ref.csv"))
write_csv(ela_sens_fn(dem_ref)$elasticity %>% as_tibble(), 
          paste0("analysis/demographic_output/",model,"/elas_ref.csv"))





#==========
#========== Variance Partitioning
#==========

# new function to process sensitivity output
part_process = function(x){
  data.frame(
    chain = x$chain,
    step = x$step,
    s_rec = x$sensitivity[1,4],
    s_s21 = x$sensitivity[2,1],
    s_s32 = x$sensitivity[3,2],
    s_s43 = x$sensitivity[4,3],
    s_s44 = x$sensitivity[4,4]
  )
}

# function to calcualte covariance matrix and partition variances
cov_fn = function(data, sensdata){
  d = data %>% 
    select(rec, s21, s32, s43, s44) %>%
    as.matrix()
  cv = cov(d)
  
  sensd = 
    sensdata %>%
    rename(
      rec = s_rec, s21 = s_s21, s32 = s_s32, s43 = s_s43, s44 = s_s44
    ) %>%
    select(rec, s21, s32, s43, s44) %>%
    as.numeric()
  
  cvmat = hadamard.prod(cv, sensd %*% t(sensd))
  rowSums(cvmat)/sum(cvmat)
}

# define mean matrix for calculating sensitivities
dem_pars_mean = dem_pars %>%
  group_by(chain, step) %>%
  summarize(
    rec = mean(rec),
    s21 = mean(s21),
    s32 = mean(s32),
    s43 = mean(s43),
    s44 = mean(s44)
  ) %>%
  mutate(
    year = 1986
  ) %>%
  ungroup()


# calculate sensitivities
sensMat = parallel::mclapply(1:nrow(dem_pars_mean), 
                FUN=function(x){
                  dem_pars_mean[x,] %>%
                    ela_sens_fn %>%
                    part_process
                }) %>%
  bind_rows %>%
  as_tibble()

# create matrix and chains and steps
var_par = dem_pars %>%
  expand(nesting(chain, step)) %>%
  mutate(row = 1:length(chain)) %>%
  split(.$row) %>%
  parallel::mclapply(function(x){
    dem_pars_ = dem_pars %>%
      filter(chain == x$chain, step == x$step)
    sensMat_ = sensMat %>%
      filter(chain == x$chain, step == x$step)
    as_tibble(t(cov_fn(dem_pars_, sensMat_)))
  }) %>%
  bind_rows() %>% 
  gather(name, val, 1:5) %>%
  group_by(name) %>%
  summarize(lower16 = quantile(val, 0.18),
            middle = median(val),
            upper84 = quantile(val, 0.84))

# export
write_csv(var_par, paste0("analysis/demographic_output/",model,"/var_par.csv"))





functions{
  // define autocorrelation function
  real acor_fn(int n, int y, matrix resid) {
    real resid_mean; // residual mean
    real acov[y-1]; // autocovariance
    real avar[y-1]; // variance
    real acor; // autocorrelation
    resid_mean = mean(resid[n,]);
    for (t in 2:y){
      acov[t-1] = (resid[n,t] - resid_mean)*(resid[n,t-1] - resid_mean); // autocovariance
      avar[t-1] = (resid[n,t-1] - resid_mean)*(resid[n,t-1] - resid_mean); // variance
      }
      acor = mean(acov[])/mean(avar[]); // autocorrelation
      return acor;
    }
}
data{
  // declare variables
  int nYears; // number of years
  int nStages; // number of stages
  matrix[nStages, nYears] N_obs; // observed abundance
  matrix[nStages, nYears] n_obs; // log observed abundance
}
transformed data{
  // declare variables
  real scale; // scaling factor for per capita recruitment (mean abundance of age 2 fish)
  real sd_obs; // sd of observed counts on log scale 
  scale = mean(N_obs[2,])/mean(N_obs[4,]);
  sd_obs = sd(n_obs[2:4,]);
}
parameters{
  real log_rho; // log_rho
  real logit_phi_init[nStages-1]; // initial value for logit_phi
  real<lower=0> sig_phi; // sd for stochastic change in log_phi
  real z_phi[nStages-1, nYears-2]; // random variate for stochastic change in log_phi
  real<lower=0> sig_obs; // sd for observation error
  real<lower=0> N_init[nStages]; // initial log abundance
}
transformed parameters{
  real<lower=0> N[nStages, nYears]; // estimated abundance
  real logit_phi[nStages-1, nYears-1]; // surival probility (on logit scale)
  real rho; // per recruitment
  real phi[nStages-1, nYears-1]; //survival probability
  real n_pred[nStages, nYears]; // estimated abundance on log scale
  // initial values
  for(i in 1:nStages){
    N[i,1] = N_init[i];
  }
  for(i in 1:(nStages-1)){
      logit_phi[i,1] = logit_phi_init[i];
    }
  // time-varying parameters
  for(t in 2:(nYears - 1)){
    for(i in 1:(nStages-1)){
      logit_phi[i,t] = logit_phi[i,t-1] + sig_phi*z_phi[i,t-1];
    }
  // tranform parameters
  rho = exp(log_rho);
  phi = inv_logit(logit_phi);
  }
  // project demographic process
  for (t in 2:nYears) {
    N[1,t] = scale*rho*N[4,t-1];
    N[2,t] = N[1,t-1];
    N[3,t] = phi[1,t-1]*N[2,t-1];
    N[4,t] = phi[2,t-1]*N[3,t-1] + phi[3,t-1]*N[4,t-1];
  }
  // log predicted values
  for (i in 1:nStages){
    n_pred[i,] = log(N[i,]);
  }
}
model{
  //random deviates for random walks
  for (i in 1:(nStages-1)){
    z_phi[i,] ~ normal(0, 1);
  }
  // priors
  log_rho ~ normal(0, 1);
  logit_phi_init ~ normal(0, 1);
  sig_phi ~ normal(0, 2) T[0,];
  sig_obs ~ normal(0, 1) T[0, ];
  for(i in 1:nStages){
    N_init[i] ~ lognormal(n_obs[i,1], sig_obs*sd_obs);
  }
  // observation process
  for(i in 2:nStages){
    for (t in 1:nYears) {
    N_obs[i,t] ~ lognormal(n_pred[i,t], sig_obs*sd_obs);
    }
  }
}
generated quantities {
  real N_tot[nYears]; // total abundance
  real rho_scale; // rho on natural scale
  real lambda[nYears-1]; // realized per capita growth rate
  real sig_obs_sd; // observation error on natural scale
  matrix[nStages-1, nYears] n_sim; // simulated log-abundance
  matrix[nStages-1, nYears] resid_obs; // observed residual
  matrix[nStages-1, nYears] resid_sim; // simulated residual 
  real acor_obs[nStages-1]; // observed residual autocorrelation
  real acor_sim[nStages-1]; // simulated residual autocorrelation
  matrix[nStages-1,nYears] log_lik; // log-likelihood
  real log_lik_sum;
  for (t in 1:nYears){
      N_tot[t] = sum(N[,t]);
      for (i in 2:nStages){
        n_sim[i-1,t] = normal_rng(n_pred[i,t], sig_obs*sd_obs);
        resid_obs[i-1,t] = n_obs[i,t] - n_pred[i,t];
        resid_sim[i-1,t] = n_sim[i-1,t] - n_pred[i,t];
        log_lik[i-1,t] = lognormal_lpdf(N_obs[i,t]|n_pred[i,t], sig_obs*sd_obs);
      }
    }
  log_lik_sum = sum(log_lik);
  rho_scale = scale*rho; 
  for (t in 1:(nYears-1)){
    lambda[t] = N_tot[t+1]/N_tot[t];
    sig_obs_sd = sig_obs*sd_obs;
  }
  for (i in 1:(nStages-1)){
    acor_obs[i] = acor_fn(i, nYears, resid_obs);
    acor_sim[i] = acor_fn(i, nYears, resid_sim);
  }
}

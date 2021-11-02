//=======================================================================================

// ARCTIC CHARR MATRIX POPULATION MODEL WITH TIME-VARYING PARAMETERS

//=======================================================================================





//=======================================================================================

data {
  
  // declare parameters
  int n_stages;                                    // number of stages
  int n_sites;                                     // number of sites
  int n_years;                                     // number of years
  int y[n_stages, n_sites, n_years];               // survey data
  real xp[7, 2];                                   // values for prior parameterization
  int<lower=0, upper=1> tvr[2];                    // binary for time-varying rates 
  int slsm[4];                                     // mapping of survival sd's
  real k;                                          // scaling for population density 
  int n_zero[n_stages, n_years];
  int n_nonzero[n_stages, n_years];
}

//=======================================================================================

parameters {
  
  // declare parameters
  real lr0;                                        // initial recruitment
  real ls0[n_stages];                              // initial survival
  real<lower=0> x0[n_stages];                      // initial population density
  real<lower=0> slr[tvr[1]];                       // sd for recruitment random walk
  real<lower=0> sls[tvr[2] * max(slsm)];                     // sd for survival random walk 
  real zr[tvr[1] * (n_years - 1)];                 // deviates for recruitment random walk
  real zs[tvr[2] * n_stages, 
          tvr[2] * (n_years - 1)];                 // deviates for survival random walk   
  real<lower=0, upper=1> p[n_stages];              // detection probability  
  real<lower=0, upper=1> theta;                    // bernoulli zero probability

}

//=======================================================================================

transformed parameters {
  
  // declare parameterrs
  real lr[tvr[1] * (n_years - 2) + 1];             // log recruitment         
  real r[tvr[1] * (n_years - 2) + 1];              // recruitment
  real ls[n_stages, tvr[2] * (n_years - 2) + 1];   // logit survival       
  real s[n_stages, tvr[2] * (n_years - 2) + 1];    // survival          
  matrix[n_stages, n_years] x;                     // population density         
  
  
  
  // initial values
  lr[1] = lr0;
  for (i in 1 : n_stages) {
    ls[i, 1] = ls0[i];
    x[i, 1] = x0[i];
  }
  
  
  
  // time-varying demographic rates
  if (tvr[1] == 1 || tvr[2] == 1) {
    for (t in 2 : (n_years - 1)) {
    
      // recruitment
      if (tvr[1] == 1) {
        lr[t] = lr[t - 1] + slr[1] * zr[t - 1];
      }
  
      // survival
      if (tvr[2] == 1) {
        for (i in 1 : n_stages) {
          ls[i, t] = ls[i, t - 1] + sls[slsm[i]] * zs[i, t - 1];
        }
      }
      
    }
  }
  
  
  
  // tranform parameters
  r = exp(lr);
  s = inv_logit(ls);
  
  
  
  // population projection
  for (t in 2 : n_years) {
    
    // declare local variables
    real rr;
    real ss[n_stages];
    
    // recruitment
    if (tvr[1] == 0) {
      rr = r[1];
    } 
    if (tvr[1] == 1) {
      rr = r[t - 1];
    }
    
    // survival
    for (i in 1 : n_stages) {
      if (tvr[2] == 0) {
        ss[i] = s[i, 1];
      } 
      if (tvr[2] == 1) {
        ss[i] = s[i, t - 1];
      }  
    }  
  
    // projection
    x[1, t] = rr * x[4, t - 1];
    x[2, t] = ss[1] * x[1, t - 1];
    x[3, t] = ss[2] * x[2, t - 1];
    x[4, t] = ss[3] * x[3, t - 1] + ss[4] * x[4, t - 1];

  } //t
  
}

//=======================================================================================

model {
  
  // initial values
  x0 ~ exponential(1 / xp[1, 1]);
  lr0 ~ normal(xp[2, 1], xp[2, 2]);
  ls0 ~ normal(xp[3, 1], xp[3, 2]);
  
  
  
  // deviates for random walks
  if (tvr[1] == 1) {
    zr ~ normal(0, 1);
  }
  if (tvr[2] == 1) {
    for(i in 1 : n_stages) {
      zs[i] ~ normal(0, 1);
    }
  }
  
  
  
  // priors
  p ~ beta(xp[4, 1], xp[4, 2]);
  slr ~ gamma(xp[5, 1], xp[5, 1] / xp[5, 2]);
  sls ~ gamma(xp[6, 1], xp[6, 1] / xp[6, 2]);
  theta ~ beta(xp[7, 1], xp[7, 2]);
  
  
  
  // likelihood
  for (t in 1 : n_years) {
    for (i in 1 : n_stages) {
      target += n_zero[i, t] * log_sum_exp(bernoulli_lpmf(1 | theta),
                                  bernoulli_lpmf(0 | theta) + 
                                    poisson_lpmf(0 | k * p[i] * x[i, t]));
      target += n_nonzero[i, t] * bernoulli_lpmf(0 | theta);
      target += poisson_lpmf(y[i, 1 : n_nonzero[i, t] ,t] | k * p[i] * x[i, t]);
    }
  }
  
}

//=======================================================================================

generated quantities {

  // declare variables
  real log_lik [n_stages * n_years];               // stage x year log-likelihood
  real log_lik_sum;                                // total log-likelihood

  // log-likelihood by stage x year
  {
    int pos = 1;
    for (t in 1 : n_years) {
      for (i in 1 : n_stages) {
        log_lik[pos] = n_zero[i, t] * log_sum_exp(bernoulli_lpmf(1 | theta),
                                  bernoulli_lpmf(0 | theta) +
                                    poisson_lpmf(0 | k * p[i] * x[i, t])) +
                                     n_nonzero[i, t] * bernoulli_lpmf(0 | theta) +
                                      poisson_lpmf(y[i, 1 : n_nonzero[i, t] ,t] | 
                                                    k * p[i] * x[i, t]);
        pos += 1;
      }
    }
  }

  // total log-likelihood
  log_lik_sum = sum(log_lik);

}

//=======================================================================================

//=======================================================================================

// Mass x age structred demography for Myvatn arctic charr

//=======================================================================================




//=======================================================================================


functions {
  
  // Growth rate as function of size class
  vector growth(real b0_, real b1_, vector q_) {
    
    // Return growth rate as negative exponential of mass
    return exp(b0_ - b1_ * q_);
  }
  
  
  
  // Fill growth transition matrix
  matrix growth_trans(vector b_) {
    
    // Declare variables
    int s_ = num_elements(b_);
    matrix [s_, s_] m_ = rep_matrix(0, s_, s_);
    
    // Loop over columns; fill subdiagonal
    for (i in 1:(s_ - 1)) {
      m_[i + 1, i] = b_[i];
    } // i
    
    // Return
    return m_;
  }
  
  
  
  // Fill mortality transition matrix
  matrix mort_trans(vector u_) {
    
    // Fill mortalities on diagonal
    return diag_matrix(u_);
  }
  


  // Expanded mass transition matrix
  matrix trans_expand(matrix m_, int a_) {
    
    // Declare variations
    int p_ = 1;
    int s_ = rows(m_);
    matrix [a_ * s_, a_ * s_] MM_ = rep_matrix(0, a_ * s_, a_ * s_);
    
    // Expand transition matrix by age
    for (i in 1:a_) {
      MM_[p_:(p_ + s_ - 1), p_:(p_ + s_ - 1)] = m_;
      p_ += s_;
    }
    
    // Return
    return MM_;
  }
  
  
  
  // Expanded fertility matrix
  matrix fert_expand(int s_, int a_, real f, int fa_, vector q_) {
    
    // Declare variables
    matrix [s_, s_] F_ = rep_matrix(0, s_, s_);
    matrix [a_ * s_, a_ * s_] FF_ = rep_matrix(0, a_ * s_, a_ * s_);
    int p_ = 1;
    
    
    // Mass-specific fertility matrix (linear function of mass)
    for (i in 1:s_) {
      F_[1, i] = f * q_[i];
    }
    
    // Exand fertility matrix by age
    for (i in 1:a_) {
      if (i >= fa_) FF_[p_:(p_ + s_ - 1), p_:(p_ + s_ - 1)] = F_;
      p_ += s_;
    }
    
    // Return
    return FF_;
    
  }  
  
  
  
  // Convert transition rate matrix to transition probability matrix
  matrix trans_prob(matrix m0_) {
    
    // Declare variables
    int k_ = cols(m0_);
    vector [k_] v_;
    matrix [k_, k_] m1_;
    matrix [k_, k_] m2_;
    matrix [k_, k_] m3_;
    matrix [k_, k_] m4_;
    
    // Original transition matrix
    m1_ = m0_; 
    
    // Replace diagonal elements of transition rate matrix with 0s
    for (i in 1:k_) {
      m1_[i, i] = 0; 
    }
    
    // Sum rates for each starting state
    for (i in 1:k_) {
      v_[i] = sum(m0_[,i]); 
    } 
    
    // Diagonal matrix with summed rates
    m2_ = diag_matrix(v_); 
    
    // Fill diagonal of transition rate matrix with negative summed rates
    m3_ = m1_ - m2_; 
    
    // Exponentiate transition matrix (solution to transition probability ODE)
    m4_ = matrix_exp(m3_); 
    
    // Return
    return m4_;
  }
  
}


//=======================================================================================


data {
  
  // Declare variables
  int times; // number of time steps
  int masses; // number of masses (mass classes)
  int ages; // number of age classes
  vector [masses] q; // stage mass covariate
  int fa; // age class at first reproduction
  matrix [ages * masses, ages * masses] K; // vec-permutation matrix
  matrix [ages * masses, ages * masses] D; // age-transition matrix
  matrix [ages * masses, ages * masses] H; // age-fertility matrix
  int y [ages * masses, times]; // observed abundance
  int yo [ages * masses]; // index of observation status
  real p [8, 2]; // priors
  real ws; // weighting for penalization
  
}


//=======================================================================================


transformed data {
  
  // Declare variables
  matrix [ages * masses, ages * masses] Kt; // transposed vec-permutation matrix
  
  // Transpose vec-permutation matrix
  Kt = K';
  
}


//=======================================================================================


parameters {
  
  // Declare variables
  // vector <lower=0> [times] b0; // baseline growth rate
  // vector <lower=0> [times] b1; // decline in growth rate with mass
  real <lower=0> b0; // baseline growth rate
  real <lower=0> b1; // decline in growth rate with mass
  vector <lower=0> [masses] u; // mass-specific mortality rate
  real <lower=0> b0s; // random walk standard deviation for b0
  real <lower=0> b1s; // random walk standard deviation for b1
  real <lower=0> us; // random walk standard deviation for u
  real <lower=0> fs; // random walk standard deviation for f
  vector <lower=0> [times - 1] f; // mass-specific fertility
  // real <lower=0> ys; // observation error sd
  vector <lower=0> [ages * masses] x0; // initial abundance 
  
}


//=======================================================================================


transformed parameters {
  
  // Declare variables
  matrix <lower = 0> [ages * masses, times] x; // abundance
  real <lower = 0> AA [ages * masses, ages * masses, times - 1]; // projection matrix array
  

  // Initial abundance
  x[, 1] = x0; // initial abundance
  
  
  { 
    
    // Declare local variables
    matrix [masses, masses] M; // mass-transition matrix
    matrix [ages * masses, ages * masses] MM; // age-expanded mass-transition matrix
    matrix [ages * masses, ages * masses] FF; // age-expanded mass-fertility matrix
    matrix [ages * masses, ages * masses] A; // projection matrix
    matrix [ages * masses, ages * masses] Mp; // age x mass transition matrix
    matrix [ages * masses, ages * masses] Fp; // age x mass fertility matrix
    
    for (t in 2:times) {
  
      // Mass-specific ransition matrix (common across all ages)
      M = trans_prob(growth_trans(growth(b0, b1, q)) + 
                                    mort_trans(u));
  
      // Expand transition matrix by age
      MM = trans_expand(M, ages);
  
      // Expand fertility matrix by age
      FF = fert_expand(masses, ages, f[t - 1], fa, q);
      
      // Vec-permutation
      Mp = Kt * D * K * MM;
      Fp = Kt * H * K * FF;
        
      // Projection matrix
      A = (Mp + Fp);
      AA[, , t - 1] = to_array_2d(A);
      
      // Project dynamics
      x[, t] = A * x[, t - 1]; 
  
    } // t
    
  }
  
}


//=======================================================================================


model {
  
  // Random walk standard deviations
  // b0s ~ gamma(p[1, 1], p[1, 2]);
  // b1s ~ gamma(p[2, 1], p[2, 2]);
  fs ~ gamma(p[2, 1], p[2, 2]);
  // us ~ gamma(p[3, 1], p[3, 2]);
  
  // Initial values for random walks
  b0 ~ gamma(p[4, 1], p[4, 2]);
  b1 ~ gamma(p[5, 1], p[5, 2]);
  u ~ gamma(p[6, 1], p[6, 2]);
  
  // Mass-specific fertility
  f[1] ~ gamma(p[7, 1], p[7, 2]);
  
  // Random walks
  for (t in 2:(times - 1)) {
    
    // Growth rates
    // b0[t] ~ normal(b0[t - 1], b0s) T[0, ];
    // b1[t] ~ normal(b1[t - 1], b0s) T[0, ];
     f[t] ~ normal(f[t - 1], fs) T[0, ];
    
    // Mortality 
    // for (i in 1:(masses - 1)) {
    //   u[i, t] ~ normal(u[i, t - 1], us) T[0, ];
    // } // i
  } // t
  
  // Observation error sd
  // ys ~ gamma(p[8, 1], p[8, 2]);
  
  // Likelihood
  for (i in 1:ages * masses) {
    if (yo[i] == 1) {
      x0[i] ~ gamma(p[8, 1], p[8, 1]/(y[i, 1] + 0.1));
      y[i, ] ~ poisson(x[i, ] + 0.1); 
    }
    if (yo[i] == 0) {
      x0[i] ~ gamma(p[8, 1], p[8, 2]);
    }
  }

  // {
  //   real w;
  //   real s;
  //   for (i in 1:(masses - 1)) {
  //     w = mean(u[i, ]);
  //     s = ws * us;
  //     target += normal_lpdf(u[i, ] | w, s) - 1/(s * sqrt(2 * pi()));
  //   } // i
  // }
  
}

//=======================================================================================



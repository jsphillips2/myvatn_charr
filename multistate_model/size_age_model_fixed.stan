//=======================================================================================

// Mass x age structred demography for Myvatn arctic charr

//=======================================================================================




//=======================================================================================


functions {
  
  // Growth rate as function of size class
  vector growth(real b0_, real b1_, vector q_) {
    
    // Return growth rate as negative exponential of mass
    return b0_ * exp(-b1_ * q_);
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
  matrix [ages * masses, times] y; // observed abundance
  int yo [ages * masses]; // index of observation status
  real p [5, 2]; // priors
  
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
  real <lower=0> b0; // baseline growth rate
  real <lower=0> b1; // decline in growth rate with mass
  vector <lower=0> [masses] u; // mass-specific mortality rate
  real <lower=0> f; // mass-specific fertility
  real <lower=0> ys; // observation error sd
  vector <lower=0> [ages * masses] x0; // initial abundance 
  
}


//=======================================================================================


transformed parameters {
  
  // Declare variables
  matrix [ages * masses, times] x; // abundance
  real AA [ages * masses, ages * masses, times - 1]; // projection matrix array
  

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
      M = trans_prob(growth_trans(growth(b0, b1, q)) + mort_trans(u));
  
      // Expand transition matrix by age
      MM = trans_expand(M, ages);
  
      // Expand fertility matrix by age
      FF = fert_expand(masses, ages, f, fa, q);
      
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
  
  // Baseline growth rate
  b0 ~ gamma(p[1, 1], p[1, 2]);
  
  // Decline in growth rate with mass
  b1 ~ gamma(p[2, 1], p[2, 2]);
  
  // Mass-specific mortality rate
  u ~ gamma(p[3, 1], p[3, 2]);
  
  // Mass-specific fertility
  f ~ gamma(p[4, 1], p[4, 2]);
  
  // Observation error sd
  ys ~ gamma(p[5, 1], p[5, 2]);
  
  // Initial abundance
  for (i in 1:ages * masses) {
    x0[i] ~ normal(y[i, 1], ys) T[0, ];
  }
  
  // Likelihood
  for (i in 1:ages * masses) {
    if (yo[i]  == 1) {
     y[i, ] ~ normal(x[i, ], ys); 
    }
  }
  
}

//=======================================================================================



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
    int s_ = num_elements(b_) + 1;
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
  int stages; // number of stages (mass classes)
  int ages; // number of age classes
  vector [stages] q; // mean stage masses
  int fa; // age class at first reproduction
  matrix [ages * stages, ages * stags] K; // vec-permutation matrix
  matrix [ages * stages, ages * stags] D; // age-transition matrix
  matrix [ages * stages, ages * stags] H; // age-fertility matrix
  
}


//=======================================================================================


transformed data {
  
  // Declare variables
  matrix [ages * stages, ages * stags] Kt; // transposed vec-permutation matrix
  
  // Transpose vec-permutation matrix
  Kt = K';
  
}


//=======================================================================================


parameters {
  
}


//=======================================================================================


transformed parameters {
  
  // Declare variables
  matrix [ages * stages, times] x; // abundance
  real [ages * stages, ages * stages, times - 1] AA; // projection matrix array
  

  // Initial abundance
  x[, 1] = x0; // initial abundance
  
  
  { 
    
    // Declare local variables
    matrix [stages, stages] M; // mass-transition matrix
    matrix [ages * stages, ages * stages] MM; // age-expanded mass-transition matrix
    matrix [ages * stages, ages * stages] FF; // age-expanded mass-fertility matrix
    matrix [ages * stages, ages * stages] A; // projection matrix
    
    for (t in 2:times) {
  
      // Mass-specific ransition matrix (common across all ages)
      M = trans_prob(growth_trans(b0[t - 1], b1[t - 1], q) + mort_trans(u[t - 1, ]));
  
      // Expand transition matrix by age
      MM = trans_expand(M, ages);
  
      // Expand fertility matrix by age
      FF = fert_expand(stages, ages, f[t - 1], fa, q);
      
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
  
}


//=======================================================================================


generated quantities {
  
}


//=======================================================================================


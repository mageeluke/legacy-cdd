data {
  int<lower=0> K;   		 		// number of fixed effects
  int<lower=0> N;                		// num data points
  matrix[N, K] x;  		 		// fixed effects predictor matrix
  int<lower=0,upper=1> y[N];     		// binary survival outcomes (0 = death, 1 = survival)
  int<lower=1> SK;               		// num species predictors
  int<lower=1> SJ;               		// num species
  int<lower=1,upper=SJ> sp[N];   		// species IDs
  matrix[N, SK] SX;              		// species predictors
  int<lower=1> nquad;          			// num quadrats
  int<lower=1,upper=nquad> quadnum[N];		// Quadrat random effect
}
parameters {
  matrix[SK, SJ] Sz;
  cholesky_factor_corr[SK] SL_Omega;
  vector<lower=0,upper=pi()/2>[SK] Stau_unif;
  vector[K] beta;	                            	// fixed effects
  real<lower=0> sigma_QUAD;    				// Quadrat error 
  vector[nquad] z_QUAD;      				// Quadrat values (as z-scores)
}
transformed parameters {
  matrix[SJ, SK] SB;
  matrix[SK, SK] sigma_SB;
  vector<lower=0>[SK] Stau;    	 
  vector[nquad] QUAD;      				// Actual quadrat IDs					
  for (k in 1:SK) Stau[k] = 2.5 * tan(Stau_unif[k]);
  SB = (diag_pre_multiply(Stau,SL_Omega) * Sz)';					// matrix of species random effects
  sigma_SB = diag_pre_multiply(Stau, SL_Omega) * diag_pre_multiply(Stau, SL_Omega)';	// calculate the covariance matrix
  QUAD = z_QUAD * sigma_QUAD;				// Calculation of actual quadrat values for quadrat random intercept
}
model {
  to_vector(Sz) ~ std_normal();
  SL_Omega ~ lkj_corr_cholesky(1);
  beta ~ normal(0, 10);
  sigma_QUAD ~ gamma(2,1);	// 1/A, A = mode of prior, A = 1
  z_QUAD ~ std_normal();
  y ~ bernoulli_logit((x * beta) + rows_dot_product(SB[sp], SX) + QUAD[quadnum]);
}
generated quantities {
  vector[N] fitted;
  vector[N] log_lik;
  for (n in 1:N) {
    fitted[n] = (x[n,] * beta) + dot_product(SB[sp[n],], SX[n,]) + QUAD[quadnum[n]];
  }
  for (n in 1:N) {
    log_lik[n] = bernoulli_logit_lpmf(y[n] | fitted[n]);
  }
}

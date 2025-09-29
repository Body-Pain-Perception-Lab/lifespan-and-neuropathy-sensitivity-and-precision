// Stan model for 2ifc simulated data.
// This model is used to elicit priors for future experiments
// PF formulation : Weibull
// Author: Arthur S. Courtin  
// License: MIT (see LICENSE file) 

data{
  int N_cdt;                //Total number of cold detection trials
  int N_wdt;                //Total number of warm detection trials
  
  int P;                    //Total number of participants
  
  row_vector[N_cdt] x_cdt;  //Stimulus intensities for each cold detection trial
  row_vector[N_wdt] x_wdt;  //Stimulus intensities for each warm detection trial
  
  array[N_cdt] int y_cdt;   //Number of detections for each cold detection intensity
  array[N_wdt] int y_wdt;   //Number of detections for each warm detection intensity
  
  array[N_cdt] int n_cdt;   //Number of presentations for each cold detection intensity
  array[N_wdt] int n_wdt;   //Number of presentations for each warm detection intensity
  
  array[N_cdt] int p_cdt;   //Participant index for each cold detection trial
  array[N_wdt] int p_wdt;   //Participant index for each warm detection trial
}
transformed data{
  int C=2;                       //Number of conditions
  int I=2;                       //Number of parameters of interest
  int U=1;                       //Number of nuisance parameters
  int K=(C*I+U);                 //Number of coefficients
}
parameters{
  vector[K] mu;
  vector<lower=0>[K] tau;
  matrix[K,P] z;
  cholesky_factor_corr[K] L;
}
transformed parameters{
  matrix[K,P] delta;              //Participant deviation from group mean
  matrix[K,P] kappa;              //Intercepts with random participant effect
  array[C] row_vector[P] alpha;   //Thresholds
  array[C] row_vector[P] beta;    //Slopes
  row_vector[P] lambda;           //Lapse rates
  row_vector[N_cdt] theta_cdt;    //Response probability for cold detection trials
  row_vector[N_wdt] theta_wdt;    //Response probability for warm detection trials

  delta = diag_pre_multiply(tau, L) * z;
  for(idx in 1:K){
    kappa[idx,] = mu[idx] + delta[idx,];
  }
  
  alpha[1] = exp(kappa[1,]);
  alpha[2] = exp(kappa[2,]);
  beta[1] = exp(kappa[3,]);
  beta[2] = exp(kappa[4,]);
  lambda = .5*inv_logit(kappa[5,]);

  theta_cdt = 1-exp(-(x_cdt./alpha[1][p_cdt]).^beta[1][p_cdt]);
  theta_wdt = 1-exp(-(x_wdt./alpha[2][p_wdt]).^beta[2][p_cdt]);
}
model{
  //Priors
  mu[1:2] ~ normal(0,1);
  mu[3:4] ~ normal(1.5,1.5);
  mu[5] ~ normal(-4,1);

  tau[1:2] ~ normal(0,1);
  tau[3:4] ~ normal(0,1.5);
  tau[5] ~ normal(0,1);

  L ~ lkj_corr_cholesky(2);
  
  to_vector(z) ~ std_normal();

  //Likelihoods  
  y_cdt ~ binomial(n_cdt,0.5 + (0.5 - lambda[p_cdt]) .* theta_cdt);
  y_wdt ~ binomial(n_wdt,0.5 + (0.5 - lambda[p_wdt]) .* theta_wdt);
}
generated quantities{
  matrix[K,K] omega = L * L'; //Correlation matrix
}

// Stan model for the ageing and neuropathy dataset but without predictors
// This model is used to elicit priors for future studies
// PF formulation : Gaussian
// Author: Arthur S. Courtin  
// License: MIT (see LICENSE file) 


data{
  int N_cdt;                //Total number of cold detection trials
  int N_wdt;                //Total number of warm detection trials
  int N_cpt;                //Total number of cold pain trials
  int N_hpt;                //Total number of heat pain trials
  
  int P;                    //Total number of participants
  
  row_vector[N_cdt] x_cdt;  //Stimulus intensities for each cold detection trial
  row_vector[N_wdt] x_wdt;  //Stimulus intensities for each warm detection trial
  row_vector[N_cpt] x_cpt;  //Stimulus intensities for each cold pain trial
  row_vector[N_hpt] x_hpt;  //Stimulus intensities for each heat pain trial
  
  array[N_cdt] int y_cdt;   //Number of detections for each cold detection intensity
  array[N_wdt] int y_wdt;   //Number of detections for each warm detection intensity
  array[N_cpt] int y_cpt;   //Number of detections for each cold pain intensity
  array[N_hpt] int y_hpt;   //Number of detections for each heat pain intensity
  
  array[N_cdt] int n_cdt;   //Number of presentations for each cold detection intensity
  array[N_wdt] int n_wdt;   //Number of presentations for each warm detection intensity
  array[N_cpt] int n_cpt;   //Number of presentations for each cold pain intensity
  array[N_hpt] int n_hpt;   //Number of presentations for each heat pain intensity
  
  array[N_cdt] int p_cdt;   //Participant index for each cold detection trial
  array[N_wdt] int p_wdt;   //Participant index for each warm detection trial
  array[N_cpt] int p_cpt;   //Participant index for each cold pain trial
  array[N_hpt] int p_hpt;   //Participant index for each heat pain trial
}
transformed data{
  int C=4;                       //Number of conditions
  int I=2;                       //Number of parameters of interest
  int U=2;                       //Number of nuisance parameters
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
  matrix[K,P] kappa;            //Intercepts with random participant effect
  array[C] row_vector[P] alpha;   //Thresholds
  array[C] row_vector[P] beta;    //Slopes
  row_vector[P] gamma;            //Guess rates
  row_vector[P] lambda;           //Lapse rates
  row_vector[N_cdt] theta_cdt;    //Response probability for cold detection trials
  row_vector[N_wdt] theta_wdt;    //Response probability for warm detection trials
  row_vector[N_cpt] theta_cpt;    //Response probability for cold pain trials
  row_vector[N_hpt] theta_hpt;    //Response probability for heat pain trials

  delta = diag_pre_multiply(tau, L) * z;
  for(idx in 1:K){
    kappa[idx,] = mu[idx] + delta[idx,];
  }
  
  alpha[1] = kappa[1,];
  alpha[2] = kappa[2,];
  alpha[3] = kappa[3,];
  alpha[4] = kappa[4,];
  beta[1] = exp(kappa[5,]);
  beta[2] = exp(kappa[6,]);
  beta[3] = exp(kappa[7,]);
  beta[4] = exp(kappa[8,]);
  gamma = .5*inv_logit(kappa[9,]);
  lambda = .5*inv_logit(kappa[10,]);

  theta_cdt = Phi(beta[1][p_cdt].*(x_cdt-alpha[1][p_cdt]));
  theta_wdt = Phi(beta[2][p_wdt].*(x_wdt-alpha[2][p_wdt]));
  theta_cpt = Phi(beta[3][p_cpt].*(x_cpt-alpha[3][p_cpt]));
  theta_hpt = Phi(beta[4][p_hpt].*(x_hpt-alpha[4][p_hpt]));
}
model{
  //Priors
  mu[1] ~ normal(0.67,0.09*4);
  mu[2] ~ normal(4.33,0.31*4);
  mu[3] ~ normal(22.48,1.18*4);
  mu[4] ~ normal(15.11,0.42*4);
  mu[5] ~ normal(1.19,0.22*4);
  mu[6] ~ normal(-0.23,0.15*4);
  mu[7] ~ normal(-1.60,0.11*4);
  mu[8] ~ normal(-0.22,0.15*4);
  mu[9] ~ normal(-3.33,0.40);
  mu[10] ~ normal(-3.84,0.36);

  tau[1] ~ normal(0.54,0.09*4);
  tau[2] ~ normal(1.90,0.24*4);
  tau[3] ~ normal(7.39,0.83*4);
  tau[4] ~ normal(2.67,0.30*4);
  tau[5] ~ normal(1.11,0.23*4);
  tau[6] ~ normal(0.73,0.13*4);
  tau[7] ~ normal(0.50,0.13*4);
  tau[8] ~ normal(0.79,0.14*4);
  tau[9] ~ normal(0.66,0.35);
  tau[10] ~ normal(1.13,0.45);

  L ~ lkj_corr_cholesky(2);
  
  to_vector(z) ~ std_normal();

  //Likelihoods  
  y_cdt ~ binomial(n_cdt,gamma[p_cdt] + (1 - gamma[p_cdt] - lambda[p_cdt]) .* theta_cdt);
  y_wdt ~ binomial(n_wdt,gamma[p_wdt] + (1 - gamma[p_wdt] - lambda[p_wdt]) .* theta_wdt);
  y_cpt ~ binomial(n_cpt,gamma[p_cpt] + (1 - gamma[p_cpt] - lambda[p_cpt]) .* theta_cpt);
  y_hpt ~ binomial(n_hpt,gamma[p_hpt] + (1 - gamma[p_hpt] - lambda[p_hpt]) .* theta_hpt);
}
generated quantities{
  matrix[K,K] omega = L * L'; //Correlation matrix
}

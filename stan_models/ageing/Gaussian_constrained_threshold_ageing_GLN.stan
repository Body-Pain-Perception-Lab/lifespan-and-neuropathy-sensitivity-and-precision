// Stan model for the ageing and neuropathy dataset
// PF formulation : Gaussian
// Age-threshold relationship formulation : log-linear
// Guess and lapse rate treatment: nuisance (random intercepts but no age effects)
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
  
  array[N_cdt] int y_cdt;   //Binary response (0 not detected/1 detected) for each cold detection trial
  array[N_wdt] int y_wdt;   //Binary response (0 not detected/1 detected) for each warm detection trial
  array[N_cpt] int y_cpt;   //Binary response (0 not burning/1 burning) for each cold pain trial
  array[N_hpt] int y_hpt;   //Binary response (0 not burning/1 burning) for each heat pain trial
  
  array[N_cdt] int p_cdt;   //Participant index for each cold detection trial
  array[N_wdt] int p_wdt;   //Participant index for each warm detection trial
  array[N_cpt] int p_cpt;   //Participant index for each cold pain trial
  array[N_hpt] int p_hpt;   //Participant index for each heat pain trial
  
  row_vector[P] a;          //Age in years of each participant
}
transformed data{
  int C=4;                       //Number of conditions
  int I=2;                       //Number of parameters of interest
  int U=2;                       //Number of nuisance parameters
  int K=(C*I+U);                 //Number of coefficients

  row_vector[P] ra=(a-min(a));   //Transformed age
}
parameters{
  vector[K+C*I] mu;
  vector<lower=0>[K] tau;
  matrix[K,P] z;
  cholesky_factor_corr[K] L;
}
transformed parameters{
  matrix[K,P] delta;              //Participant deviation from group mean
  matrix[K,P] kappa_0;            //Intercepts with random participant effect
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
    kappa_0[idx,] = mu[idx] + delta[idx,];
  }
  
  alpha[1] = exp(kappa_0[1,]+mu[K+1].*ra);
  alpha[2] = exp(kappa_0[2,]+mu[K+2].*ra);
  alpha[3] = exp(kappa_0[3,]+mu[K+3].*ra);
  alpha[4] = exp(kappa_0[4,]+mu[K+4].*ra);
  beta[1] = exp(kappa_0[5,]+mu[K+5].*ra);
  beta[2] = exp(kappa_0[6,]+mu[K+6].*ra);
  beta[3] = exp(kappa_0[7,]+mu[K+7].*ra);
  beta[4] = exp(kappa_0[8,]+mu[K+8].*ra);
  gamma = .5*inv_logit(kappa_0[9,]);
  lambda = .5*inv_logit(kappa_0[10,]);

  theta_cdt = Phi(beta[1][p_cdt].*(x_cdt-alpha[1][p_cdt]));
  theta_wdt = Phi(beta[2][p_wdt].*(x_wdt-alpha[2][p_wdt]));
  theta_cpt = Phi(beta[3][p_cpt].*(x_cpt-alpha[3][p_cpt]));
  theta_hpt = Phi(beta[4][p_hpt].*(x_hpt-alpha[4][p_hpt]));
}
model{
  //Priors
  //Kappa_0
  mu[1] ~ normal(-0.52,0.15*4);
  mu[2] ~ normal(1.35,0.06*4);
  mu[3] ~ normal(3.06,0.08*4);
  mu[4] ~ normal(2.71,0.03*4);
  mu[5] ~ normal(1.20,0.20*4);
  mu[6] ~ normal(-0.24,0.15*4);
  mu[7] ~ normal(-1.64,0.12*4);
  mu[8] ~ normal(-0.30,0.16*4);
  mu[9] ~ normal(-3.32,0.41);
  mu[10] ~ normal(-3.74,0.35);
  
  //Kappa_1
  mu[11] ~ normal(0,0.92/60);
  mu[12] ~ normal(0,0.53/60);
  mu[13] ~ normal(0,0.38/60);
  mu[14] ~ normal(0,0.17/60);
  mu[15] ~ normal(0,1.09/60);
  mu[16] ~ normal(0,0.75/60);
  mu[17] ~ normal(0,0.53/60);
  mu[18] ~ normal(0,0.85/60);
  
  //Kappa_0
  tau[1] ~ normal(0.92,0.13*4);
  tau[2] ~ normal(0.53,0.07*4);
  tau[3] ~ normal(0.38,0.04*4);
  tau[4] ~ normal(0.17,0.04*4);
  tau[5] ~ normal(1.09,0.19*4);
  tau[6] ~ normal(0.75,0.14*4);
  tau[7] ~ normal(0.53,0.11*4);
  tau[8] ~ normal(0.85,0.14*4);
  tau[9] ~ normal(0.63,0.36);
  tau[10] ~ normal(0.55,0.38);

  L ~ lkj_corr_cholesky(2);
  
  to_vector(z) ~ std_normal();

  //Likelihoods  
  y_cdt ~ bernoulli(gamma[p_cdt] + (1 - gamma[p_cdt] - lambda[p_cdt]) .* theta_cdt);
  y_wdt ~ bernoulli(gamma[p_wdt] + (1 - gamma[p_wdt] - lambda[p_wdt]) .* theta_wdt);
  y_cpt ~ bernoulli(gamma[p_cpt] + (1 - gamma[p_cpt] - lambda[p_cpt]) .* theta_cpt);
  y_hpt ~ bernoulli(gamma[p_hpt] + (1 - gamma[p_hpt] - lambda[p_hpt]) .* theta_hpt);
}
generated quantities{
  matrix[K,K] omega = L * L'; //Correlation matrix

  vector[N_cdt] log_lik_cdt; //Log-likelihood for each cold detection trial
  vector[N_wdt] log_lik_wdt; //Log-likelihood for each warm detection trial
  vector[N_cpt] log_lik_cpt; //Log-likelihood for each cold pain trial
  vector[N_hpt] log_lik_hpt; //Log-likelihood for each heat pain trial
  vector[N_cdt+N_wdt+N_cpt+N_hpt] log_lik; //Log-likelihood for each trial

  for(n in 1:N_cdt){
    log_lik_cdt[n] = bernoulli_lpmf(y_cdt[n]|gamma[p_cdt[n]] + (1 - gamma[p_cdt[n]] - lambda[p_cdt[n]]) .* theta_cdt[n]);
  }
  for(n in 1:N_wdt){
    log_lik_wdt[n] = bernoulli_lpmf(y_wdt[n]|gamma[p_wdt[n]] + (1 - gamma[p_wdt[n]] - lambda[p_wdt[n]]) .* theta_wdt[n]);
  }
  for(n in 1:N_cpt){
    log_lik_cpt[n] = bernoulli_lpmf(y_cpt[n]|gamma[p_cpt[n]] + (1 - gamma[p_cpt[n]] - lambda[p_cpt[n]]) .* theta_cpt[n]);
  }
  for(n in 1:N_hpt){
    log_lik_hpt[n] = bernoulli_lpmf(y_hpt[n]|gamma[p_hpt[n]] + (1 - gamma[p_hpt[n]] - lambda[p_hpt[n]]) .* theta_hpt[n]);
  }
  vector[N_cdt + N_wdt]  temp1 = append_row(log_lik_cdt, log_lik_wdt);
  vector[N_cpt + N_hpt]  temp2 = append_row(log_lik_cpt, log_lik_hpt);
  log_lik = append_row(temp1, temp2);
}

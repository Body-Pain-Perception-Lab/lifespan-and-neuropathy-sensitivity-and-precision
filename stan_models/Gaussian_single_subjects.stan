// Stan model for the discrimination analyes using the ageing and neuropathy dataset
// PF formulation : Gaussian
// This model does not take explanatory variables, it only estimate PF values at the subject level.
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
}
transformed data{
  int C=4;                       //Number of conditions
}
parameters{
  array[C] row_vector[P] alpha;       //Thresholds
  array[C] row_vector[P] log_beta;    //log-Slopes
  array[C] row_vector[P] logit_gamma; //logit-transformed guess rates
  array[C] row_vector[P] logit_lambda;//logit-transformed lapse rates
}
transformed parameters{
  array[C] row_vector[P] beta;    //Slopes
  array[C] row_vector[P] gamma;   //Guess rates
  array[C] row_vector[P] lambda;  //Lapse rates
  
  row_vector[N_cdt] theta_cdt;    //Response probability for cold detection trials
  row_vector[N_wdt] theta_wdt;    //Response probability for warm detection trials
  row_vector[N_cpt] theta_cpt;    //Response probability for cold pain trials
  row_vector[N_hpt] theta_hpt;    //Response probability for heat pain trials

  beta = exp(log_beta);
  gamma = inv_logit(logit_gamma);
  lambda = inv_logit(logit_lambda);

  theta_cdt = .5 * gamma[1][p_cdt] + (1 - .5 * gamma[1][p_cdt] - .5 * lambda[1][p_cdt]) .* Phi(beta[1][p_cdt].*(x_cdt-alpha[1][p_cdt]));
  theta_wdt = .5 * gamma[2][p_wdt] + (1 - .5 * gamma[2][p_wdt] - .5 * lambda[2][p_wdt]) .* Phi(beta[2][p_wdt].*(x_wdt-alpha[2][p_wdt]));
  theta_cpt = .5 * gamma[3][p_cpt] + (1 - .5 * gamma[3][p_cpt] - .5 * lambda[3][p_cpt]) .* Phi(beta[3][p_cpt].*(x_cpt-alpha[3][p_cpt]));
  theta_hpt = .5 * gamma[4][p_hpt] + (1 - .5 * gamma[4][p_hpt] - .5 * lambda[4][p_hpt]) .* Phi(beta[4][p_hpt].*(x_hpt-alpha[4][p_hpt]));
}
model{
  //Priors
  //eta_0_0
  alpha[1] ~ normal(0.67,sqrt(0.09^2+0.54^2));
  alpha[2] ~ normal(4.33,sqrt(0.31^2+1.90^2));
  alpha[3] ~ normal(22.48,sqrt(1.18^2+7.39^2));
  alpha[4] ~ normal(15.11,sqrt(0.42^2+2.67^2));
  
  log_beta[1] ~ normal(1.19,sqrt(0.22^2+1.11^2));
  log_beta[2] ~ normal(-0.23,sqrt(0.15^2+0.73^2));
  log_beta[3] ~ normal(-1.60,sqrt(0.11^2+0.50^2));
  log_beta[4] ~ normal(-0.22,sqrt(0.15^2+0.79^2));
  
  logit_gamma[1] ~ normal(-3.33,0.40);
  logit_gamma[2] ~ normal(-3.33,0.40);
  logit_gamma[3] ~ normal(-3.33,0.40);
  logit_gamma[4] ~ normal(-3.33,0.40);
  
  logit_lambda[1] ~ normal(-3.84,0.36);
  logit_lambda[2] ~ normal(-3.84,0.36);
  logit_lambda[3] ~ normal(-3.84,0.36);
  logit_lambda[4] ~ normal(-3.84,0.36);
  
  //Likelihoods  
  y_cdt ~ bernoulli(theta_cdt);
  y_wdt ~ bernoulli(theta_wdt);
  y_cpt ~ bernoulli(theta_cpt);
  y_hpt ~ bernoulli(theta_hpt);
}


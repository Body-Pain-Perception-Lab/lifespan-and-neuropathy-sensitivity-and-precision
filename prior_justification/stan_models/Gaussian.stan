data{
  int N_cdt;
  int N_wdt;
  int N_cpt;
  int N_hpt;
  int P;
  
  row_vector[N_cdt] x_cdt;
  row_vector[N_wdt] x_wdt;
  row_vector[N_cpt] x_cpt;
  row_vector[N_hpt] x_hpt;
  
  array[N_cdt] int y_cdt;
  array[N_wdt] int y_wdt;
  array[N_cpt] int y_cpt;
  array[N_hpt] int y_hpt;
  
  array[N_cdt] int p_cdt;
  array[N_wdt] int p_wdt;
  array[N_cpt] int p_cpt;
  array[N_hpt] int p_hpt;
}
transformed data{
  int C=4;
  int I=2;
  int U=2;
  int K=C*I+U;
}
parameters{
  vector[K] mu;
  vector<lower=0>[K] tau;
  matrix[K,P] z;
  cholesky_factor_corr[K] L;
}
transformed parameters{
  matrix[K,P] delta;
  matrix[K,P] kappa;
  array[C] row_vector[P] alpha;
  array[C] row_vector[P] beta;
  row_vector[P] gamma;
  row_vector[P] lambda;
  row_vector[N_cdt] theta_cdt;
  row_vector[N_wdt] theta_wdt;
  row_vector[N_cpt] theta_cpt;
  row_vector[N_hpt] theta_hpt;

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

  theta_cdt =  Phi(beta[1][p_cdt] .* (x_cdt - alpha[1][p_cdt]));
  theta_wdt =  Phi(beta[2][p_wdt] .* (x_wdt - alpha[2][p_wdt]));
  theta_cpt =  Phi(beta[3][p_cpt] .* (x_cpt - alpha[3][p_cpt]));
  theta_hpt =  Phi(beta[4][p_hpt] .* (x_hpt - alpha[4][p_hpt]));
}
model{
  mu[1:2] ~ normal(5,5);
  mu[3:4] ~ normal(15,10);
  mu[5:8] ~ normal(0,1);
  mu[9:10] ~ normal(-4,0.5);

  tau[1:2] ~ normal(0,5);
  tau[3:4] ~ normal(0,10);
  tau[5:8] ~ normal(0,1);
  tau[9:10] ~ normal(0,0.5);
  
  L ~ lkj_corr_cholesky(2);
  
  to_vector(z) ~ std_normal();
  
  y_cdt ~ bernoulli(gamma[p_cdt] + (1 - gamma[p_cdt] - lambda[p_cdt]) .* theta_cdt);
  y_wdt ~ bernoulli(gamma[p_wdt] + (1 - gamma[p_wdt] - lambda[p_wdt]) .* theta_wdt);
  y_cpt ~ bernoulli(gamma[p_cpt] + (1 - gamma[p_cpt] - lambda[p_cpt]) .* theta_cpt);
  y_hpt ~ bernoulli(gamma[p_hpt] + (1 - gamma[p_hpt] - lambda[p_hpt]) .* theta_hpt);
}
generated quantities{
  vector[N_cdt] log_lik_cdt;
  vector[N_wdt] log_lik_wdt;
  vector[N_cpt] log_lik_cpt;
  vector[N_hpt] log_lik_hpt;
  vector[N_cdt+N_wdt+N_cpt+N_hpt] log_lik;

  for(n in 1:N_cdt){
    log_lik_cdt[n] = bernoulli_logit_lpmf(y_cdt[n]|gamma[p_cdt[n]] + (1 - gamma[p_cdt[n]] - lambda[p_cdt[n]]) .* theta_cdt[n]);
  }
  for(n in 1:N_wdt){
    log_lik_wdt[n] = bernoulli_logit_lpmf(y_wdt[n]|gamma[p_wdt[n]] + (1 - gamma[p_wdt[n]] - lambda[p_wdt[n]]) .* theta_wdt[n]);
  }
  for(n in 1:N_cpt){
    log_lik_cpt[n] = bernoulli_logit_lpmf(y_cpt[n]|gamma[p_cpt[n]] + (1 - gamma[p_cpt[n]] - lambda[p_cpt[n]]) .* theta_cpt[n]);
  }  
  for(n in 1:N_hpt){
    log_lik_hpt[n] = bernoulli_logit_lpmf(y_hpt[n]|gamma[p_hpt[n]] + (1 - gamma[p_hpt[n]] - lambda[p_hpt[n]]) .* theta_hpt[n]);
  }
  log_lik[1:N_cdt] = log_lik_cdt;
  log_lik[N_cdt+1:N_cdt+N_wdt] = log_lik_wdt;
  log_lik[N_cdt+N_wdt+1:N_cdt+N_wdt+N_cpt] = log_lik_cpt;
  log_lik[N_cdt+N_wdt+N_cpt+1:N_cdt+N_wdt+N_cpt+N_hpt] = log_lik_hpt;
}

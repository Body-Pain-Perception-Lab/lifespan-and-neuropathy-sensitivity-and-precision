// Stan model for the ageing and neuropathy dataset
// PF formulation : Gaussian
// Age-threshold relationship formulation : log-linear
// Guess and lapse rate treatment: full parametrization
// Age - neuropathy intercation : yes
// Random slope for neuropathy : yes
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
  row_vector[P] h;          //Status of each participant (0 controls/1 patients)
}
transformed data{
  int C=4;                       //Number of conditions
  int I=2;                       //Number of parameters of interest
  int U=2;                       //Number of nuisance parameters
  int H=2;                       //Number of status
  int K=(C*I+U)*2*H;             //Number of main effects
  int T=(C*I+U);                 //Number of random effects

  row_vector[P] ra=(a-min(a));   //Transformed age
}
parameters{
  vector[K] mu;              //Group mean effects
  vector<lower=0>[T] tau;    //Between participant SD
  matrix[T,P] z;             //Participant random effects
  cholesky_factor_corr[T] L; //Cholesky factor of the correlation matrix
}
transformed parameters{
  matrix[T,P] delta;              //Participant deviation from group mean
  matrix[T,P] eta_0_0;        //All participants intercept
  array[C] row_vector[P] alpha;   //Thresholds
  array[C] row_vector[P] beta;    //Slopes
  row_vector[P] gamma;            //Guess rates
  row_vector[P] lambda;           //Lapse rates
  row_vector[N_cdt] theta_cdt;    //Response probability for cold detection trials
  row_vector[N_wdt] theta_wdt;    //Response probability for warm detection trials
  row_vector[N_cpt] theta_cpt;    //Response probability for cold pain trials
  row_vector[N_hpt] theta_hpt;    //Response probability for heat pain trials

  delta = diag_pre_multiply(tau, L) * z;
  for(idx in 1:T){
    eta_0_0[idx,] = mu[idx] + delta[idx,];
  }
  
  alpha[1] = exp(eta_0_0[1,] + mu[T+1] .* h + (mu[2*T+1] + mu[3*T+1] .* h) .*ra);
  alpha[2] = exp(eta_0_0[2,] + mu[T+2] .* h + (mu[2*T+2] + mu[3*T+2] .* h) .*ra);
  alpha[3] = exp(eta_0_0[3,] + mu[T+3] .* h + (mu[2*T+3] + mu[3*T+3] .* h) .*ra);
  alpha[4] = exp(eta_0_0[4,] + mu[T+4] .* h + (mu[2*T+4] + mu[3*T+4] .* h) .*ra);
  beta[1] = exp(eta_0_0[5,] + mu[T+5] .* h + (mu[2*T+5] + mu[3*T+5] .* h) .*ra);
  beta[2] = exp(eta_0_0[6,] + mu[T+6] .* h + (mu[2*T+6] + mu[3*T+6] .* h) .*ra);
  beta[3] = exp(eta_0_0[7,] + mu[T+7] .* h + (mu[2*T+7] + mu[3*T+7] .* h) .*ra);
  beta[4] = exp(eta_0_0[8,] + mu[T+8] .* h + (mu[2*T+8] + mu[3*T+8] .* h) .*ra);
  gamma = .5*inv_logit(eta_0_0[9,] + mu[T+9] .* h + (mu[2*T+9] + mu[3*T+9] .* h) .*ra);
  lambda = .5*inv_logit(eta_0_0[10,] + mu[T+10] .* h + (mu[2*T+10] + mu[3*T+10] .* h) .*ra);

  theta_cdt = Phi(beta[1][p_cdt].*(x_cdt-alpha[1][p_cdt]));
  theta_wdt = Phi(beta[2][p_wdt].*(x_wdt-alpha[2][p_wdt]));
  theta_cpt = Phi(beta[3][p_cpt].*(x_cpt-alpha[3][p_cpt]));
  theta_hpt = Phi(beta[4][p_hpt].*(x_hpt-alpha[4][p_hpt]));
}
model{
  //Priors
  //Eta_0_0
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
  //Eta_0_1
  mu[11] ~ normal(0,0.92);
  mu[12] ~ normal(0,0.53);
  mu[13] ~ normal(0,0.38);
  mu[14] ~ normal(0,0.17);
  mu[15] ~ normal(0,1.09);
  mu[16] ~ normal(0,0.75);
  mu[17] ~ normal(0,0.53);
  mu[18] ~ normal(0,0.85);
  mu[19] ~ normal(0,0.63);
  mu[20] ~ normal(0,0.55);
  //Eta_1_0
  mu[21] ~ normal(0,0.92/60);
  mu[22] ~ normal(0,0.53/60);
  mu[23] ~ normal(0,0.38/60);
  mu[24] ~ normal(0,0.17/60);
  mu[25] ~ normal(0,1.09/60);
  mu[26] ~ normal(0,0.75/60);
  mu[27] ~ normal(0,0.53/60);
  mu[28] ~ normal(0,0.85/60);
  mu[29] ~ normal(0,0.63/60);
  mu[30] ~ normal(0,0.55/60);
  //Eta_1_1
  mu[31] ~ normal(0,0.92/60);
  mu[32] ~ normal(0,0.53/60);
  mu[33] ~ normal(0,0.38/60);
  mu[34] ~ normal(0,0.17/60);
  mu[35] ~ normal(0,1.09/60);
  mu[36] ~ normal(0,0.75/60);
  mu[37] ~ normal(0,0.53/60);
  mu[38] ~ normal(0,0.85/60);
  mu[39] ~ normal(0,0.63/60);
  mu[40] ~ normal(0,0.55/60);
  
  //Eta_0_0
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

  L ~ lkj_corr_cholesky(1);
  
  to_vector(z) ~ std_normal();
  
  //Likelihoods  
  y_cdt ~ bernoulli(gamma[p_cdt] + (1 - gamma[p_cdt] - lambda[p_cdt]) .* theta_cdt);
  y_wdt ~ bernoulli(gamma[p_wdt] + (1 - gamma[p_wdt] - lambda[p_wdt]) .* theta_wdt);
  y_cpt ~ bernoulli(gamma[p_cpt] + (1 - gamma[p_cpt] - lambda[p_cpt]) .* theta_cpt);
  y_hpt ~ bernoulli(gamma[p_hpt] + (1 - gamma[p_hpt] - lambda[p_hpt]) .* theta_hpt);
}
generated quantities{
  matrix[T,T] omega = L * L'; //Correlation matrix
}

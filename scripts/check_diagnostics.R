# Script print the different sampling diagnostics for the different models
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Ageing ####
model_names<-c(
  "fits/ageing/Gaussian_unconstrained_threshold_ageing_GLN.rds",
  "fits/ageing/Gaussian_constrained_threshold_ageing_GLN.rds",
  "fits/ageing/Quick_constrained_threshold_ageing_GLN.rds"
)
for(m in 1:3){
  fit<-readRDS(model_names[m])
  print(model_names[m])
  print(fit$diagnostic_summary())
  
  s<-fit$summary(variables = c('mu'))
  print('mu')
  print(max(s$rhat))
  print(min(s$ess_bulk))
  print(min(s$ess_bulk))
  
  s<-fit$summary(variables = c('tau'))
  print('tau')
  print(max(s$rhat))
  print(min(s$ess_bulk))
  print(min(s$ess_bulk))
  
  s<-fit$summary(variables = c('alpha','beta'))
  print('alpha and beta')
  print(max(s$rhat))
  print(min(s$ess_bulk))
  print(min(s$ess_bulk))
  
  s<-fit$summary(variables = c('lambda','gamma'))
  print('gamma and lambda')
  print(max(s$rhat))
  print(min(s$ess_bulk))
  print(min(s$ess_bulk))
}

#### Ageing and neuropathy ####
model_names<-c(
  "fits/ageing_and_neuropathy/Gaussian_unconstrained_threshold_ageing_and_neuropathy_GLN_NI_H.rds",
  "fits/ageing_and_neuropathy/Gaussian_constrained_threshold_ageing_and_neuropathy_GLN_NI_H.rds",
  "fits/ageing_and_neuropathy/Quick_constrained_threshold_ageing_and_neuropathy_GLN_NI_H.rds"
)
for(m in 1:3){
  fit<-readRDS(model_names[m])
  print(model_names[m])
  print(fit$diagnostic_summary())
  
  s<-fit$summary(variables = c('mu'))
  print('mu')
  print(max(s$rhat))
  print(min(s$ess_bulk))
  print(min(s$ess_bulk))
  
  s<-fit$summary(variables = c('tau'))
  print('tau')
  print(max(s$rhat))
  print(min(s$ess_bulk))
  print(min(s$ess_bulk))
  
  s<-fit$summary(variables = c('alpha','beta'))
  print('alpha and beta')
  print(max(s$rhat))
  print(min(s$ess_bulk))
  print(min(s$ess_bulk))
  
  s<-fit$summary(variables = c('lambda','gamma'))
  print('gamma and lambda')
  print(max(s$rhat))
  print(min(s$ess_bulk))
  print(min(s$ess_bulk))
}

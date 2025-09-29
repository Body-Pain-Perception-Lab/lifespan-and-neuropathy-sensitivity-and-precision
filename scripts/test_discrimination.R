# Script to test discrimination accuracy for different scores including various PF parameters.
# Two version of the analysis exist, one that uses the posterior means for the different hyperparameters
# and one that fully propagates uncertainty to the discrimination analysis.
# The later is reported in the manuscript 
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set up environement ####
library(tidyverse)
library(bayesplot)
library(mvtnorm)

#### Retrieve fit ####
fit_an<-readRDS("results/fits/ageing_and_neuropathy/Gaussian_unconstrained_threshold_ageing_and_neuropathy_GLN_NI_H.rds")
draws<-fit_an$draws(variables = c('mu','tau','omega'),format='df')


#### Posterior mode ####
set.seed(12345)
M<-1000 #number of simulations
N<-1000 #sample size in each group for each simulation

#Design matrix for the different scores
DM<-matrix(
  c(
    1,rep(0,0),1,rep(0,7),
    1,rep(0,1),1,rep(0,6),
    1,rep(0,0),1,1,rep(0,6),
    
    1,rep(0,2),1,rep(0,5),
    1,rep(0,3),1,rep(0,4),
    1,rep(0,2),1,1,rep(0,4),
    
    1,rep(0,4),1,rep(0,3),
    1,rep(0,5),1,rep(0,2),
    1,rep(0,4),1,1,rep(0,2),
    
    1,rep(0,6),1,rep(0,1),
    1,rep(0,7),1,rep(0,0),
    1,rep(0,6),1,1,rep(0,0),
    
    1,rep(1,8), rep(0,0)
  ),
  ncol = 9,
  byrow = T
)
D<-dim(DM)[1]

G<-100

status<-rep(0:1,N)

#Create empty data frames and fill them with simulations
s_and_s<-tibble(
  sensitivity=rep(NA,M*D*(G+1)),
  specificity=rep(NA,M*D*(G+1)),
  model=rep(NA,M*D*(G+1)),
  simulation=rep(NA,M*D*(G+1))
)
AUC<-tibble(
  AUC=rep(NA,M*D),
  model=rep(NA,M*D),
  simulation=rep(NA,M*D)
)
for(mdx in 1:M){
  print(mdx)
  
  #sample train set group parameters
  d<-sample(1:8000,1)
  mu_train<-numeric(26)
  for(idx in 1:length(mu_train)){
    mu_train[idx]<-draws %>% pull(paste0("mu[", idx, "]")) %>% mean()
  }
  tau_train<-numeric(10)
  for(idx in 1:length(tau_train)){
    tau_train[idx]<-draws %>% pull(paste0("tau[", idx, "]")) %>% mean()
  }
  omega_train<-matrix(,10,10)
  for(idx in 1:length(tau_train)){
    for(jdx in 1:length(tau_train)){
      omega_train[idx,jdx]<-draws %>% pull(paste0("omega[", idx,',',jdx, "]")) %>% mean()
      omega_train[jdx,idx]<-draws %>% pull(paste0("omega[", idx,',',jdx, "]")) %>% mean()
    }
  }
  sigma_train<-diag(tau_train)%*%omega_train%*%diag(tau_train)
  
  #sample train set individual parameters
  p_train<-rmvnorm(2*N,mu_train[1:length(tau_train)],sigma_train)
  age_train<-runif(2*N,20,60)
  
  parameters_train<-tibble(
    age=rep(NA,2*8*N),
    status=rep(NA,2*8*N),
    parameter=rep(NA,2*8*N),
    participant=rep(NA,2*8*N),
    value=rep(NA,2*8*N)
  )
  for(parameter in 1:8){
    parameters_train$value[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-p_train[,parameter]+mu_train[10+parameter]*status+mu_train[18+parameter]*age_train
    parameters_train$age[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-20+age_train
    parameters_train$status[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-status
    parameters_train$parameter[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-parameter
    parameters_train$participant[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-sort(rep(1:N,2))
  }
  
  #sample test set group parameters
  d<-sample(1:8000,1)
  mu_test<-numeric(26)
  for(idx in 1:length(mu_train)){
    mu_test[idx]<-draws %>% pull(paste0("mu[", idx, "]")) %>% mean()
  }
  tau_test<-numeric(10)
  for(idx in 1:length(tau_train)){
    tau_test[idx]<-draws %>% pull(paste0("tau[", idx, "]")) %>% mean()
  }
  omega_test<-matrix(,10,10)
  for(idx in 1:length(tau_train)){
    for(jdx in 1:length(tau_train)){
      omega_test[idx,jdx]<-draws %>% pull(paste0("omega[", idx,',',jdx, "]")) %>% mean()
      omega_test[jdx,idx]<-draws %>% pull(paste0("omega[", idx,',',jdx, "]")) %>% mean()
    }
  }
  sigma_test<-diag(tau_test)%*%omega_test%*%diag(tau_test)
  
  #sample test set individual parameters
  p_test<-rmvnorm(2*N,mu_test[1:length(tau_test)],sigma_test)
  age_test<-runif(2*N,20,60)
  
  parameters_test<-tibble(
    age=rep(NA,2*8*N),
    status=rep(NA,2*8*N),
    parameter=rep(NA,2*8*N),
    participant=rep(NA,2*8*N),
    value=rep(NA,2*8*N)
  )
  for(parameter in 1:8){
    parameters_test$value[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-p_test[,parameter]+mu_test[10+parameter]*status+mu_test[18+parameter]*age_test
    parameters_test$age[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-20+age_test
    parameters_test$status[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-status
    parameters_test$parameter[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-parameter
    parameters_test$participant[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-sort(rep(1:N,2))
  }
  
  #assess discrimination for each model
  for(model in 1:D){ 
    #neutralize unwanted parameters
    data_train<-
      parameters_train %>% 
      mutate(parameter=factor(parameter,labels=c('cdt','wdt','cpt','hpt','cds','wds','cps','hps'))) %>% 
      pivot_wider(names_from = parameter,values_from = value) %>% 
      mutate(
        cdt=DM[model,2]*cdt,
        cds=DM[model,3]*cds,
        wdt=DM[model,4]*wdt,
        wds=DM[model,5]*wds,
        cpt=DM[model,6]*cpt,
        cps=DM[model,7]*cps,
        hpt=DM[model,8]*hpt,
        hps=DM[model,9]*hps,
        age=age
      )
    #fit logistic regression on the train set to obtain coefficients that will be used to compute scores
    logistic_model<-glm(
      status~
        cdt+wdt+cpt+hpt+cds+wds+cps+hps+
        cdt:age+wdt:age+cpt:age+hpt:age+cds:age+wds:age+cps:age+hps:age,
      data = data_train,
      family = 'binomial'
    )
    for (idx in 1:length(logistic_model$coefficients)){
      if(is.na(logistic_model$coefficients[idx])){
        logistic_model$coefficients[idx]<-0
      }
    }
    coeffs<-logistic_model$coefficients
    
    #compute scores for the test set
    scores<-
      parameters_test %>% 
      mutate(parameter=factor(parameter,labels=c('cdt','wdt','cpt','hpt','cds','wds','cps','hps'))) %>% 
      pivot_wider(names_from = parameter,values_from = value) %>% 
      mutate(
        cdt=DM[model,2]*cdt,
        cds=DM[model,3]*cds,
        wdt=DM[model,4]*wdt,
        wds=DM[model,5]*wds,
        cpt=DM[model,6]*cpt,
        cps=DM[model,7]*cps,
        hpt=DM[model,8]*hpt,
        hps=DM[model,9]*hps,
        age=age
      )%>% 
      mutate(score=
               coeffs[1]+
               coeffs[2]*cdt+
               coeffs[3]*wdt+
               coeffs[4]*cpt+
               coeffs[5]*hpt+
               coeffs[6]*cds+
               coeffs[7]*wds+
               coeffs[8]*cps+
               coeffs[9]*hps+
               coeffs[10]*cdt*age+
               coeffs[11]*wdt*age+
               coeffs[12]*cpt*age+
               coeffs[13]*hpt*age+
               coeffs[14]*cds*age+
               coeffs[15]*wds*age+
               coeffs[16]*cps*age+
               coeffs[17]*hps*age
      ) %>% 
      select(participant,status,score)
    
    #assess sensitivity, specificity and AUC for the scores at different cut-offs
    cut_offs<-seq(min(scores$score),max(scores$score),(max(scores$score)-min(scores$score))/G)
    sensitivity<-specificity<-vector(,G+1)
    for(idx in 1:length(cut_offs)){
      TP<-sum((scores$score>=cut_offs[idx])&(scores$status==1))
      TN<-sum((scores$score<cut_offs[idx])&(scores$status==0))
      FP<-sum((scores$score>=cut_offs[idx])&(scores$status==0))
      FN<-sum((scores$score<cut_offs[idx])&(scores$status==1))
      sensitivity[idx]<-TP/(TP+FN)
      specificity[idx]<-TN/(TN+FP)
    }
    
    s_and_s$sensitivity[(1+(model-1)*(G+1)+(mdx-1)*D*(G+1)):(G+1+(model-1)*(G+1)+(mdx-1)*D*(G+1))]<-sensitivity
    s_and_s$specificity[(1+(model-1)*(G+1)+(mdx-1)*D*(G+1)):(G+1+(model-1)*(G+1)+(mdx-1)*D*(G+1))]<-specificity
    s_and_s$model[(1+(model-1)*(G+1)+(mdx-1)*D*(G+1)):(G+1+(model-1)*(G+1)+(mdx-1)*D*(G+1))]<-model
    s_and_s$simulation[(1+(model-1)*(G+1)+(mdx-1)*D*(G+1)):(G+1+(model-1)*(G+1)+(mdx-1)*D*(G+1))]<-mdx
    
    auroc=0
    for(idx in 1:G){
      auroc=
        auroc+
        0.5*
        (sensitivity[idx]+sensitivity[idx+1])*
        (specificity[idx+1]-specificity[idx])
    }
    
    AUC$AUC[model+(mdx-1)*D]<-auroc
    AUC$model[model+(mdx-1)*D]<-model
    AUC$simulation[model+(mdx-1)*D]<-mdx
  }
}
saveRDS(s_and_s,'results/S_and_S_PM.rds')

AUC<-
  AUC %>% 
  mutate(
    model=factor(model,labels = c('CDT','CDS','CDT&CDS','WDT','WDS','WDT&WDS','CPT','CPS','CPT&CPS','HPT','HPS','HPT&HPS','All')),
  )
AUC$task[AUC$model%in%c('CDT','CDS','CDT&CDS')]<-'CD'
AUC$task[AUC$model%in%c('WDT','WDS','WDT&WDS')]<-'WD'
AUC$task[AUC$model%in%c('CPT','CPS','CPT&CPS')]<-'CP'
AUC$task[AUC$model%in%c('HPT','HPS','HPT&HPS')]<-'HP'
AUC$param[AUC$model%in%c('CDT','WDT','CPT','HPT')]<-'Threshold'
AUC$param[AUC$model%in%c('CDS','WDS','CPS','HPS')]<-'Slope'
AUC$param[AUC$model%in%c('CDT&CDS','WDT&WDS','CPT&CPS','HPT&HPS')]<-'Both'

AUC$param<-factor(AUC$param,levels=c('Threshold','Slope','Both','All'))
AUC$task<-factor(AUC$task,levels=c('CD','WD','CP','HP'))

saveRDS(AUC,'results/AUC_PM.rds')

#### Full uncertainty propagation ####
set.seed(12345)
M<-1000 #number of simulations
N<-1000 #sample size in each group for each simulation

#Design matrix for the different scores
DM<-matrix(
  c(
    1,rep(0,0),1,rep(0,7),
    1,rep(0,1),1,rep(0,6),
    1,rep(0,0),1,1,rep(0,6),
    
    1,rep(0,2),1,rep(0,5),
    1,rep(0,3),1,rep(0,4),
    1,rep(0,2),1,1,rep(0,4),
    
    1,rep(0,4),1,rep(0,3),
    1,rep(0,5),1,rep(0,2),
    1,rep(0,4),1,1,rep(0,2),
    
    1,rep(0,6),1,rep(0,1),
    1,rep(0,7),1,rep(0,0),
    1,rep(0,6),1,1,rep(0,0),
    
    1,rep(1,3), rep(0,5)
  ),
  ncol = 9,
  byrow = T
)
D<-dim(DM)[1]

G<-100

status<-rep(0:1,N)

#Create empty data frames and fill them with simulations
s_and_s<-tibble(
  sensitivity=rep(NA,M*D*(G+1)),
  specificity=rep(NA,M*D*(G+1)),
  model=rep(NA,M*D*(G+1)),
  simulation=rep(NA,M*D*(G+1))
)
AUC<-tibble(
  AUC=rep(NA,M*D),
  model=rep(NA,M*D),
  simulation=rep(NA,M*D)
)
for(mdx in 1:M){
  print(mdx)
  
  #sample train set group parameters
  d<-sample(1:8000,1)
  mu_train<-numeric(26)
  for(idx in 1:length(mu_train)){
    temp <- draws %>% pull(paste0("mu[", idx, "]"))
    mu_train[idx]<-temp[d]
  }
  tau_train<-numeric(10)
  for(idx in 1:length(tau_train)){
    temp <- draws %>% pull(paste0("tau[", idx, "]"))
    tau_train[idx]<-temp[d]
  }
  omega_train<-matrix(,10,10)
  for(idx in 1:length(tau_train)){
    for(jdx in 1:length(tau_train)){
      temp <- draws %>% pull(paste0("omega[", idx,',',jdx, "]"))
      omega_train[idx,jdx]<-temp[d]
      omega_train[jdx,idx]<-temp[d]
    }
  }
  sigma_train<-diag(tau_train)%*%omega_train%*%diag(tau_train)
  
  #sample train set individual parameters
  p_train<-rmvnorm(2*N,mu_train[1:length(tau_train)],sigma_train)
  age_train<-runif(2*N,20,60)
  
  parameters_train<-tibble(
    age=rep(NA,2*8*N),
    status=rep(NA,2*8*N),
    parameter=rep(NA,2*8*N),
    participant=rep(NA,2*8*N),
    value=rep(NA,2*8*N)
  )
  for(parameter in 1:8){
    parameters_train$value[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-p_train[,parameter]+mu_train[10+parameter]*status+mu_train[18+parameter]*age_train
    parameters_train$age[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-20+age_train
    parameters_train$status[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-status
    parameters_train$parameter[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-parameter
    parameters_train$participant[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-sort(rep(1:N,2))
  }
  
  #sample test set group parameters
  d<-sample(1:8000,1)
  mu_test<-numeric(26)
  for(idx in 1:length(mu_test)){
    temp <- draws %>% pull(paste0("mu[", idx, "]"))
    mu_test[idx]<-temp[d]
  }
  tau_test<-numeric(10)
  for(idx in 1:length(tau_test)){
    temp <- draws %>% pull(paste0("tau[", idx, "]"))
    tau_test[idx]<-temp[d]
  }
  omega_test<-matrix(,10,10)
  for(idx in 1:length(tau_test)){
    for(jdx in 1:length(tau_test)){
      temp <- draws %>% pull(paste0("omega[", idx,',',jdx, "]"))
      omega_test[idx,jdx]<-temp[d]
      omega_test[jdx,idx]<-temp[d]
    }
  }
  sigma_test<-diag(tau_test)%*%omega_test%*%diag(tau_test)
  
  #sample test set individual parameters
  p_test<-rmvnorm(2*N,mu_test[1:length(tau_test)],sigma_test)
  age_test<-runif(2*N,20,60)
  
  parameters_test<-tibble(
    age=rep(NA,2*8*N),
    status=rep(NA,2*8*N),
    parameter=rep(NA,2*8*N),
    participant=rep(NA,2*8*N),
    value=rep(NA,2*8*N)
  )
  for(parameter in 1:8){
    parameters_test$value[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-p_test[,parameter]+mu_test[10+parameter]*status+mu_test[18+parameter]*age_test
    parameters_test$age[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-20+age_test
    parameters_test$status[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-status
    parameters_test$parameter[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-parameter
    parameters_test$participant[(1+2*N*(parameter-1)):(2*N+2*N*(parameter-1))]<-sort(rep(1:N,2))
  }
  
  #assess discrimination for each model
  for(model in 1:D){ 
    #neutralize unwanted parameters
    data_train<-
      parameters_train %>% 
      mutate(parameter=factor(parameter,labels=c('cdt','wdt','cpt','hpt','cds','wds','cps','hps'))) %>% 
      pivot_wider(names_from = parameter,values_from = value) %>% 
      mutate(
        cdt=DM[model,2]*cdt,
        cds=DM[model,3]*cds,
        wdt=DM[model,4]*wdt,
        wds=DM[model,5]*wds,
        cpt=DM[model,6]*cpt,
        cps=DM[model,7]*cps,
        hpt=DM[model,8]*hpt,
        hps=DM[model,9]*hps,
        age=age
      )
    #fit logistic regression on the train set to obtain coefficients that will be used to compute scores
    logistic_model<-glm(
      status~
        cdt+wdt+cpt+hpt+cds+wds+cps+hps+
        cdt:age+wdt:age+cpt:age+hpt:age+cds:age+wds:age+cps:age+hps:age,
      data = data_train,
      family = 'binomial'
    )
    for (idx in 1:length(logistic_model$coefficients)){
      if(is.na(logistic_model$coefficients[idx])){
        logistic_model$coefficients[idx]<-0
      }
    }
    coeffs<-logistic_model$coefficients
    
    #compute scores for the test set
    scores<-
      parameters_test %>% 
      mutate(parameter=factor(parameter,labels=c('cdt','wdt','cpt','hpt','cds','wds','cps','hps'))) %>% 
      pivot_wider(names_from = parameter,values_from = value) %>% 
      mutate(
        cdt=DM[model,2]*cdt,
        cds=DM[model,3]*cds,
        wdt=DM[model,4]*wdt,
        wds=DM[model,5]*wds,
        cpt=DM[model,6]*cpt,
        cps=DM[model,7]*cps,
        hpt=DM[model,8]*hpt,
        hps=DM[model,9]*hps,
        age=age
      )%>% 
      mutate(score=
               coeffs[1]+
               coeffs[2]*cdt+
               coeffs[3]*wdt+
               coeffs[4]*cpt+
               coeffs[5]*hpt+
               coeffs[6]*cds+
               coeffs[7]*wds+
               coeffs[8]*cps+
               coeffs[9]*hps+
               coeffs[10]*cdt*age+
               coeffs[11]*wdt*age+
               coeffs[12]*cpt*age+
               coeffs[13]*hpt*age+
               coeffs[14]*cds*age+
               coeffs[15]*wds*age+
               coeffs[16]*cps*age+
               coeffs[17]*hps*age
      ) %>% 
      select(participant,status,score)
    
    #assess sensitivity, specificity and AUC for the scores at different cut-offs
    cut_offs<-seq(min(scores$score),max(scores$score),(max(scores$score)-min(scores$score))/G)
    sensitivity<-specificity<-vector(,G+1)
    for(idx in 1:length(cut_offs)){
      TP<-sum((scores$score>=cut_offs[idx])&(scores$status==1))
      TN<-sum((scores$score<cut_offs[idx])&(scores$status==0))
      FP<-sum((scores$score>=cut_offs[idx])&(scores$status==0))
      FN<-sum((scores$score<cut_offs[idx])&(scores$status==1))
      sensitivity[idx]<-TP/(TP+FN)
      specificity[idx]<-TN/(TN+FP)
    }
    
    s_and_s$sensitivity[(1+(model-1)*(G+1)+(mdx-1)*D*(G+1)):(G+1+(model-1)*(G+1)+(mdx-1)*D*(G+1))]<-sensitivity
    s_and_s$specificity[(1+(model-1)*(G+1)+(mdx-1)*D*(G+1)):(G+1+(model-1)*(G+1)+(mdx-1)*D*(G+1))]<-specificity
    s_and_s$model[(1+(model-1)*(G+1)+(mdx-1)*D*(G+1)):(G+1+(model-1)*(G+1)+(mdx-1)*D*(G+1))]<-model
    s_and_s$simulation[(1+(model-1)*(G+1)+(mdx-1)*D*(G+1)):(G+1+(model-1)*(G+1)+(mdx-1)*D*(G+1))]<-mdx
    
    auroc=0
    for(idx in 1:G){
      auroc=
        auroc+
        0.5*
        (sensitivity[idx]+sensitivity[idx+1])*
        (specificity[idx+1]-specificity[idx])
    }
    
    AUC$AUC[model+(mdx-1)*D]<-auroc
    AUC$model[model+(mdx-1)*D]<-model
    AUC$simulation[model+(mdx-1)*D]<-mdx
  }
}
saveRDS(s_and_s,'results/S_and_S_FUP.rds')

AUC<-
  AUC %>% 
  mutate(
    model=factor(model,labels = c('CDT','CDS','CDT&CDS','WDT','WDS','WDT&WDS','CPT','CPS','CPT&CPS','HPT','HPS','HPT&HPS','All')),
  )
AUC$task[AUC$model%in%c('CDT','CDS','CDT&CDS')]<-'CD'
AUC$task[AUC$model%in%c('WDT','WDS','WDT&WDS')]<-'WD'
AUC$task[AUC$model%in%c('CPT','CPS','CPT&CPS')]<-'CP'
AUC$task[AUC$model%in%c('HPT','HPS','HPT&HPS')]<-'HP'
AUC$param[AUC$model%in%c('CDT','WDT','CPT','HPT')]<-'Threshold'
AUC$param[AUC$model%in%c('CDS','WDS','CPS','HPS')]<-'Slope'
AUC$param[AUC$model%in%c('CDT&CDS','WDT&WDS','CPT&CPS','HPT&HPS')]<-'Both'

AUC$param<-factor(AUC$param,levels=c('Threshold','Slope','Both','All'))
AUC$task<-factor(AUC$task,levels=c('CD','WD','CP','HP'))

saveRDS(AUC,'results/AUC_FUP.rds')

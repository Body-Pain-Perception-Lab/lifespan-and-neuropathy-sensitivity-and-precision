# Script to fit the data simulated for the parameter and model recovery analysis
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set up environment ####
library(tidyverse)
library(cmdstanr)

if (!dir.exists("parameter_recovery/fits")) {
  dir.create("parameter_recovery/fits", recursive = TRUE)
}
if (!dir.exists("parameter_recovery/loo")) {
  dir.create("parameter_recovery/loo", recursive = TRUE)
}

model_symbols<-
  c(
    'GU',
    'GC',
    'QC'
  )
dataset_symbols<-
  c(
    'Guc',
    'Gc',
    'Qc'
  )

models=c(
  "stan_models/ageing_and_neuropathy/Gaussian_constrained_threshold_ageing_and_neuropathy_GLN_NI_H.stan",
  "stan_models/ageing_and_neuropathy/Gaussian_unconstrained_threshold_ageing_and_neuropathy_GLN_NI_H.stan",
  "stan_models/ageing_and_neuropathy/Quick_constrained_threshold_ageing_and_neuropathy_GLN_NI_H.stan"
)

#### Create functions to generate inits #####
init_f_GU<-function(){
  init<-list(
    mu=c(
      rnorm(1,0.67,0.09*4),
      rnorm(1,4.33,0.31*4),
      rnorm(1,22.48,1.18*4),
      rnorm(1,15.11,0.42*4),
      rnorm(1,1.19,0.22*4),
      rnorm(1,-0.23,0.15*4),
      rnorm(1,-1.60,0.11*4),
      rnorm(1,-0.22,0.15*4),
      rnorm(1,-3.33,0.40),
      rnorm(1,-3.84,0.36),
      
      rnorm(1,0,0.54),
      rnorm(1,0,1.90),
      rnorm(1,0,7.39),
      rnorm(1,0,2.67),
      rnorm(1,0,1.11),
      rnorm(1,0,0.73),
      rnorm(1,0,0.50),
      rnorm(1,0,0.79),
      
      rnorm(1,0,0.54/60),
      rnorm(1,0,1.90/60),
      rnorm(1,0,7.39/60),
      rnorm(1,0,2.67/60),
      rnorm(1,0,1.11/60),
      rnorm(1,0,0.73/60),
      rnorm(1,0,0.50/60),
      rnorm(1,0,0.79/60)
    ),
    tau=c(
      abs(rnorm(1,0.54,0.09*4)),
      abs(rnorm(1,1.90,0.24*4)),
      abs(rnorm(1,7.39,0.83*4)),
      abs(rnorm(1,2.67,0.30*4)),
      abs(rnorm(1,1.11,0.23*4)),
      abs(rnorm(1,0.73,0.13*4)),
      abs(rnorm(1,0.50,0.13*4)),
      abs(rnorm(1,0.79,0.14*4)),
      abs(rnorm(1,0.66,0.35)),
      abs(rnorm(1,1.13,0.45))
    )
  )
  return(init)
}
init_f_GC<-function(){
  init<-list(
    mu=c(
      rnorm(1,-0.52,0.15*4),
      rnorm(1,1.35,0.06*4),
      rnorm(1,3.06,0.08*4),
      rnorm(1,2.71,0.03*4),
      rnorm(1,1.20,0.20*4),
      rnorm(1,-0.24,0.15*4),
      rnorm(1,-1.64,0.12*4),
      rnorm(1,-0.30,0.16*4),
      rnorm(1,-3.32,0.41),
      rnorm(1,-3.74,0.35),
      
      rnorm(1,0,0.92),
      rnorm(1,0,0.53),
      rnorm(1,0,0.38),
      rnorm(1,0,0.17),
      rnorm(1,0,1.09),
      rnorm(1,0,0.75),
      rnorm(1,0,0.53),
      rnorm(1,0,0.85),
      
      rnorm(1,0,0.92/60),
      rnorm(1,0,0.53/60),
      rnorm(1,0,0.38/60),
      rnorm(1,0,0.17/60),
      rnorm(1,0,1.09/60),
      rnorm(1,0,0.75/60),
      rnorm(1,0,0.53/60),
      rnorm(1,0,0.85/60)
    ),
    tau=c(
      abs(rnorm(1,0.92,0.13*4)),
      abs(rnorm(1,0.53,0.07*4)),
      abs(rnorm(1,0.38,0.04*4)),
      abs(rnorm(1,0.17,0.04*4)),
      abs(rnorm(1,1.09,0.19*4)),
      abs(rnorm(1,0.75,0.14*4)),
      abs(rnorm(1,0.53,0.11*4)),
      abs(rnorm(1,0.85,0.14*4)),
      abs(rnorm(1,0.63,0.36)),
      abs(rnorm(1,0.55,0.38))
    )
  )
  return(init)
}
init_f_QC<-function(){
  init<-list(
    mu=c(
      rnorm(1,-0.57,0.17*4),
      rnorm(1,1.35,0.09*4),
      rnorm(1,3.09,0.06*4),
      rnorm(1,2.72,0.03*4),
      rnorm(1,0.71,0.14*4),
      rnorm(1,1.27,0.15*4),
      rnorm(1,1.51,0.11*4),
      rnorm(1,2.59,0.17*4),
      rnorm(1,-3.25,0.39),
      rnorm(1,-4.02,0.36),
      
      rnorm(1,0,1.02),
      rnorm(1,0,0.53),
      rnorm(1,0,0.41),
      rnorm(1,0,0.16),
      rnorm(1,0,0.39),
      rnorm(1,0,0.73),
      rnorm(1,0,0.33),
      rnorm(1,0,0.97),
      
      rnorm(1,0,1.02/60),
      rnorm(1,0,0.53/60),
      rnorm(1,0,0.41/60),
      rnorm(1,0,0.16/60),
      rnorm(1,0,0.39/60),
      rnorm(1,0,0.73/60),
      rnorm(1,0,0.33/60),
      rnorm(1,0,0.97/60)
    ),
    tau=c(
      abs(rnorm(1,1.02,0.15*4)),
      abs(rnorm(1,0.53,0.07*4)),
      abs(rnorm(1,0.41,0.05*4)),
      abs(rnorm(1,0.16,0.04*4)),
      abs(rnorm(1,0.39,0.19*4)),
      abs(rnorm(1,0.73,0.14*4)),
      abs(rnorm(1,0.33,0.14*4)),
      abs(rnorm(1,0.97,0.15*4)),
      abs(rnorm(1,0.86,0.36)),
      abs(rnorm(1,0.39,0.28))
    )
  )
  return(init)
}

#### Fit the simulated data ####
#loop through generative model
for (gen_idx in 1:3){
  all_trials <- read.csv(paste0("parameter_recovery/simulated_runs/",dataset_symbols[gen_idx],"_trial_info.csv"))
  all_participants <- read.csv(paste0("parameter_recovery/simulated_runs/",dataset_symbols[gen_idx],"_participant_info.csv"))
 
  #loop through simulated datasets
  for(dataset_idx in 1:15){
    print(dataset_idx)
    
    #reformat data
    participants<-
      all_participants %>%
      filter(dataset==dataset_idx,type==1)
  
    trials<-
      all_trials %>%
      filter(dataset==dataset_idx,participant %in% participants$participant)
  
    x_cdt<-trials$x[trials$type==1]
    x_wdt<-trials$x[trials$type==2]
    x_cpt<-trials$x[trials$type==3]
    x_hpt<-trials$x[trials$type==4]
  
    y_cdt<-trials$response[trials$type==1]
    y_wdt<-trials$response[trials$type==2]
    y_cpt<-trials$response[trials$type==3]
    y_hpt<-trials$response[trials$type==4]
  
    p_cdt<-trials$participant[trials$type==1]
    p_wdt<-trials$participant[trials$type==2]
    p_cpt<-trials$participant[trials$type==3]
    p_hpt<-trials$participant[trials$type==4]
  
    N_cdt<-length(x_cdt)
    N_wdt<-length(x_wdt)
    N_cpt<-length(x_cpt)
    N_hpt<-length(x_hpt)
  
    P<-length(unique(trials$participant))
  
    a<-participants$age
    h<-participants$status
  
    data<-
      list(
        x_cdt = x_cdt,
        x_wdt = x_wdt,
        x_cpt = x_cpt,
        x_hpt = x_hpt,
  
        y_cdt = y_cdt,
        y_wdt = y_wdt,
        y_cpt = y_cpt,
        y_hpt = y_hpt,
  
        p_cdt = p_cdt,
        p_wdt = p_wdt,
        p_cpt = p_cpt,
        p_hpt = p_hpt,
  
        N_cdt = N_cdt,
        N_wdt = N_wdt,
        N_cpt = N_cpt,
        N_hpt = N_hpt,
  
        P = P,
        a = a,
        h = h
      )
    #fit each model 
    for(mod_idx in 1:3){
      mod<-cmdstan_model(stan_file=models[mod_idx],
                         stanc_options = list("O1"),
                         force_recompile=T
                         )
      fit<-mod$sample(data=data,
                      seed=123,
                      chains=4,
                      parallel_chains = 4,
                      iter_warmup =2000,
                      iter_sampling= 2000,
                      max_treedepth = 12,
                      adapt_delta = .99,
                      refresh = 100,
                      init= get(paste0('init_f_', model_symbols[mod_idx]))
                      )
      fit$save_object(paste0('parameter_recovery/fits/',model_symbols[gen_idx],'_',model_symbols[mod_idx],'_',dataset_idx,'.rds'))
      
      loo<-fit$loo(cores=4,variables='log_lik')
      saveRDS(loo,paste0("parameter_recovery/loo/", model_symbols[gen_idx],'_',model_symbols[mod_idx],'_',dataset_idx,'.rds'))
      rm(loo)
      gc()
    }
  }
}



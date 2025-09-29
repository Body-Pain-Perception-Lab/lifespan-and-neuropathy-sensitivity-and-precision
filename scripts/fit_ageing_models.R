# Script to fit the different models the ageing only dataset.
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(cmdstanr)
library(tidyverse)
library(furrr)
library(loo)

Sys.setenv(PATH = "/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin")

rm(list=ls())

directory=getwd()

if (!dir.exists("sampling/inits")) {
  dir.create("sampling/inits", recursive = TRUE)
}
if (!dir.exists("sampling/c++_models")) {
  dir.create("sampling/c++_models", recursive = TRUE)
}
if (!dir.exists("sampling/output")) {
  dir.create("sampling/output", recursive = TRUE)
}
if (!dir.exists("sampling/samples")) {
  dir.create("sampling/samples", recursive = TRUE)
}
if (!dir.exists("results/fits/ageing")) {
  dir.create("results/fits/ageing", recursive = TRUE)
}
if (!dir.exists("results/loo/ageing")) {
  dir.create("results/loo/ageing", recursive = TRUE)
}

#### Define function ####
fit_model<-function(parameters){
  mod<-cmdstan_model(stan_file=paste(parameters$directory,"/stan_models/ageing/",parameters$model,sep=""),
                     dir=paste(parameters$directory,"/sampling/c++_models",sep=""),
                     stanc_options = list("O1"),
                     force_recompile=T,
                     compile_model_methods=T
  )
  
  file_name_out <- paste(parameters$directory,"/sampling/output/output_", str_remove(parameters$model,'.stan'), ".txt",sep="")
  file_name_warn <- paste(parameters$directory,"/sampling/output/output_", str_remove(parameters$model,'.stan'), "_warn.txt",sep="")
  
  con_out <- file(file_name_out, open = "wt")   # Open for writing
  con_warn <- file(file_name_warn, open = "wt") # Open for writing
  
  sink(con_out, append = TRUE)
  sink(con_warn, append = TRUE, type = "message")
  
  on.exit(sink())   
  
  print(parameters$model)
  
  mt=12
  ad=.99
  iw=2000
  is=2000
  
  fit<-mod$sample(data=parameters$data,
                  output_dir=paste(parameters$directory,"/sampling/samples",sep=""),
                  seed=12345,
                  chains=4,
                  parallel_chains = 4,
                  iter_warmup =iw,
                  iter_sampling= is,
                  save_warmup=F,
                  max_treedepth = mt,
                  adapt_delta = ad,
                  refresh = is/10,
                  init=c(paste0(parameters$directory,"/sampling/inits/init_", str_remove_all(parameters$model,".stan"), "_chain_",1:4, ".json")),
                  show_messages=T,
                  show_exceptions=T)
  
  fit$save_object(paste0(parameters$directory,"/results/fits/ageing/", str_remove_all(parameters$model,".stan"), ".rds"))

  loo<-fit$loo(cores=4,variables='log_lik')
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  loo<-fit$loo(cores=4,variables='log_lik',moment_match=T)
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/MM_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()

  loo<-fit$loo(cores=4,variables='log_lik_cdt')
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/cdt_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  loo<-fit$loo(cores=4,variables='log_lik_cdt',moment_match=T)
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/MM_cdt_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  loo<-fit$loo(cores=4,variables='log_lik_wdt')
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/wdt_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  loo<-fit$loo(cores=4,variables='log_lik_wdt',moment_match=T)
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/MM_wdt_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  loo<-fit$loo(cores=4,variables='log_lik_cpt')
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/cpt_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  loo<-fit$loo(cores=4,variables='log_lik_cpt',moment_match=T)
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/MM_cpt_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  loo<-fit$loo(cores=4,variables='log_lik_hpt')
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/hpt_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  loo<-fit$loo(cores=4,variables='log_lik_hpt',moment_match=T)
  saveRDS(loo,paste0(parameters$directory,"/results/loo/ageing/MM_hpt_", str_remove_all(parameters$model,".stan"), ".rds"))
  rm(loo)
  gc()
  sink(NULL)
  sink(NULL)  
}

possfit_model = possibly(.f = fit_model, otherwise = "Error")

init_f_Gaussian_unconstrained_threshold_ageing<-function(){
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
      
      rnorm(1,0,0.54/60),
      rnorm(1,0,1.90/60),
      rnorm(1,0,7.39/60),
      rnorm(1,0,2.67/60),
      rnorm(1,0,1.11/60),
      rnorm(1,0,0.73/60),
      rnorm(1,0,0.50/60),
      rnorm(1,0,0.79/60),
      rnorm(1,0,0.66/60),
      rnorm(1,0,1.13/60)
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
init_f_Gaussian_constrained_threshold_ageing<-function(){
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
      
      rnorm(1,0,0.92/60),
      rnorm(1,0,0.53/60),
      rnorm(1,0,0.38/60),
      rnorm(1,0,0.17/60),
      rnorm(1,0,1.09/60),
      rnorm(1,0,0.75/60),
      rnorm(1,0,0.53/60),
      rnorm(1,0,0.85/60),
      rnorm(1,0,0.63/60),
      rnorm(1,0,0.55/60)
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
init_f_Quick_constrained_threshold_ageing<-function(){
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
      
      rnorm(1,0,1.02/60),
      rnorm(1,0,0.53/60),
      rnorm(1,0,0.41/60),
      rnorm(1,0,0.16/60),
      rnorm(1,0,0.39/60),
      rnorm(1,0,0.73/60),
      rnorm(1,0,0.33/60),
      rnorm(1,0,0.97/60),
      rnorm(1,0,0.86/60),
      rnorm(1,0,0.39/60)
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

init_f_Gaussian_unconstrained_threshold_ageing_GLN<-function(){
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
init_f_Gaussian_constrained_threshold_ageing_GLN<-function(){
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
init_f_Quick_constrained_threshold_ageing_GLN<-function(){
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

#### Import and prepare data ####
d <- read_csv("data/aggregated_results.csv") %>% 
  filter(
    recording_deviates_from_target==0,
    recording_deviates_from_mean==0,
    !is.na(response),
    !subject%in%c(3019), #exclusion based on saturation check plots
    subject<5000
  ) %>% 
  arrange(subject)
participant<-unique(d$subject)

participant_info<-read_csv("data/demographics.csv") %>% 
  filter(ID%in%participant) %>% 
  arrange(ID)

P<-length(participant)
a<-participant_info$age

for(pdx in 1:P){
  d$subject[d$subject==participant[pdx]]<-pdx
}

d_cdt<-
  d %>% 
  filter(condition==1)
N_cdt<-length(d_cdt$condition)
x_cdt<-d_cdt$recorded_intensity
y_cdt<-d_cdt$response
p_cdt<-d_cdt$subject

d_wdt<-
  d %>% 
  filter(condition==2)
N_wdt<-length(d_wdt$condition)
x_wdt<-d_wdt$recorded_intensity
y_wdt<-d_wdt$response
p_wdt<-d_wdt$subject

d_cpt<-
  d %>% 
  filter(condition==3)
N_cpt<-length(d_cpt$condition)
x_cpt<-d_cpt$recorded_intensity
y_cpt<-d_cpt$response
p_cpt<-d_cpt$subject

d_hpt<-
  d %>% 
  filter(condition==4)
N_hpt<-length(d_hpt$condition)
x_hpt<-d_hpt$recorded_intensity
y_hpt<-d_hpt$response
p_hpt<-d_hpt$subject

data<-list(
  P=P,
  a=a,
  N_cdt=N_cdt,
  N_wdt=N_wdt,
  N_cpt=N_cpt,
  N_hpt=N_hpt,
  x_cdt=x_cdt,
  x_wdt=x_wdt,
  x_cpt=x_cpt,
  x_hpt=x_hpt,
  y_cdt=y_cdt,
  y_wdt=y_wdt,
  y_cpt=y_cpt,
  y_hpt=y_hpt,
  p_cdt=p_cdt,
  p_wdt=p_wdt,
  p_cpt=p_cpt,
  p_hpt=p_hpt
)

#### Prepare lists of models #############
models=c(
  'Gaussian_unconstrained_threshold_ageing.stan',
  'Gaussian_constrained_threshold_ageing.stan',
  'Quick_constrained_threshold_ageing.stan',

  'Gaussian_unconstrained_threshold_ageing_GLN.stan',
  'Gaussian_constrained_threshold_ageing_GLN.stan',
  'Quick_constrained_threshold_ageing_GLN.stan'
)

#### Prepare inits #############

for (m in models) {
  init_function_name <- paste0("init_f_", str_remove_all(m,".stan"))
  init_function <- get(init_function_name)
  
  for (c in 1:4) {
    set.seed(c)
    
    init_list <- init_function()
    
    # Define output file name
    json_file <- paste0(directory,"/sampling/inits/init_", str_remove_all(m,".stan"), "_chain_", c, ".json")
    
    # Save to JSON in a Stan-readable format
    write_stan_json(init_list, json_file)
  }
}

#### Prepare lists for fitting runs ##############
parameters_fit<-list()
for(m in 1:length(models)){
  parameters_fit[[m]]<-
    list(
      model=models[m],
      data=data,
      directory=directory
    )
}

#### Fit models (in parallel) ####
plan(multisession, workers = 3)
results<-future_map(
  parameters_fit,
  ~possfit_model(.x),
  .options = furrr_options(
    seed = T,
    scheduling = F,
    stdout = F,
    conditions = character()),
  .progress = F)











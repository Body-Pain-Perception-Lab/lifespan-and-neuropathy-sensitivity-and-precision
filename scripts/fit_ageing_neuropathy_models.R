# Script to fit the ageing model
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

if (!dir.exists("sampling/c++_models")) {
  dir.create("sampling/c++_models", recursive = TRUE)
}
if (!dir.exists("sampling/output")) {
  dir.create("sampling/output", recursive = TRUE)
}
if (!dir.exists("results/fits")) {
  dir.create("results/fits", recursive = TRUE)
}

#### Define function ####
fit_model<-function(parameters){
  file_name_out <- paste(parameters$directory,"/sampling/output/output_", str_remove(parameters$model,'.stan'), ".txt",sep="")

  con_out <- file(file_name_out, open = "wt")   # Open for writing

  sink(con_out, append = TRUE)

  on.exit(sink())   
  
  seed=12345
  mt=12
  ad=.95
  iw=2000
  is=2000  
  
  print(parameters$model)
  sprintf('Max. treedepth: %i',mt) %>% print()
  sprintf('Adapt. delta: %f',ad) %>% print()
  sprintf('Warm-up: %i',iw) %>% print()
  sprintf('Sampling: %i',is) %>% print()
  
  mod<-cmdstan_model(stan_file=paste(parameters$directory,"/stan_models/",parameters$model,sep=""),
                     dir=paste(parameters$directory,"/sampling/c++_models",sep=""),
                     stanc_options = list("O1")
  )

  pathfinder_fit<-
    mod$pathfinder(
      data=data,
      psis_resample = F,
      calculate_lp = F,
      seed = seed,
      refresh=500
    )
  
  fit<-mod$sample(data=parameters$data,
                  seed=seed,
                  chains=4,
                  parallel_chains = 4,
                  iter_warmup =iw,
                  iter_sampling= is,
                  save_warmup=F,
                  max_treedepth = mt,
                  adapt_delta = ad,
                  refresh = is/10,
                  init=pathfinder_fit,
                  show_messages=T,
                  show_exceptions=T)
  
  fit$diagnostic_summary() %>% print()
  fit$summary(c('mu','tau')) %>% print(n=40)
  
  fit$save_object(paste0(parameters$directory,"/results/fits/", str_remove_all(parameters$model,".stan"), ".rds"))
  
  sink(NULL)
}

possfit_model = possibly(.f = fit_model, otherwise = "Error")

#### Import and prepare data ####
d <- read_csv("data/aggregated_results.csv") %>% 
  filter(
    recording_deviates_from_target==0,
    recording_deviates_from_mean==0,
    !is.na(response),
    !subject%in%c(3019), #exclusion based on saturation check plots
  ) %>% 
  arrange(subject)
participant<-unique(d$subject)

participant_info<-read_csv("data/demographics.csv") %>% 
  filter(ID%in%participant) %>% 
  arrange(ID)

P<-length(participant)
a<-participant_info$age
h<-participant_info$ID>=5000

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
  h=as.numeric(h),
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
  'Gaussian_constrained_threshold_ageing_and_neuropathy_NRS_NI.stan'
  )

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
plan(multisession, workers = 1)
results<-future_map(
  parameters_fit,
  ~possfit_model(.x),
  .options = furrr_options(
    seed = T,
    scheduling = F,
    stdout = F,
    conditions = character()),
  .progress = F)











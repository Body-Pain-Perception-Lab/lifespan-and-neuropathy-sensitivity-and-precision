# Script to fit a model without predictors to the ageing and neuropathy dataset.
# The results of this fit will be used to inform priors of future yes no psi experiments.
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
if (!dir.exists("sampling/samples")) {
  dir.create("sampling/samples", recursive = TRUE)
}
if (!dir.exists("results/fits/no_predictor")) {
  dir.create("results/fits/no_predictor", recursive = TRUE)
}
#### Define function ####
init_f_Gaussian_no_predictor<-function(){
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
      rnorm(1,-3.84,0.36)
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

#### Import and prepare data ####
d <- read_csv("data/aggregated_results.csv") %>% 
  filter(
    recording_deviates_from_target==0,
    recording_deviates_from_mean==0,
    !is.na(response),
    !subject%in%c(3019) #exclusion based on saturation check plots
  ) %>% 
  arrange(subject) %>% 
  mutate(recorded_intensity=round(recorded_intensity,1)) %>% 
  group_by(subject,condition,recorded_intensity) %>% 
  summarise(y=sum(response),n=sum(!is.na(response)))
  
participant<-unique(d$subject)

P<-length(participant)

for(pdx in 1:P){
  d$subject[d$subject==participant[pdx]]<-pdx
}

d_cdt<-
  d %>% 
  filter(condition==1)
N_cdt<-length(d_cdt$condition)
x_cdt<-d_cdt$recorded_intensity
y_cdt<-d_cdt$y
n_cdt<-d_cdt$n
p_cdt<-d_cdt$subject

d_wdt<-
  d %>% 
  filter(condition==2)
N_wdt<-length(d_wdt$condition)
x_wdt<-d_wdt$recorded_intensity
y_wdt<-d_wdt$y
n_wdt<-d_wdt$n
p_wdt<-d_wdt$subject

d_cpt<-
  d %>% 
  filter(condition==3)
N_cpt<-length(d_cpt$condition)
x_cpt<-d_cpt$recorded_intensity
y_cpt<-d_cpt$y
n_cpt<-d_cpt$n
p_cpt<-d_cpt$subject

d_hpt<-
  d %>% 
  filter(condition==4)
N_hpt<-length(d_hpt$condition)
x_hpt<-d_hpt$recorded_intensity
y_hpt<-d_hpt$y
n_hpt<-d_hpt$n
p_hpt<-d_hpt$subject

data<-list(
  P=P,
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
  n_cdt=n_cdt,
  n_wdt=n_wdt,
  n_cpt=n_cpt,
  n_hpt=n_hpt,
  p_cdt=p_cdt,
  p_wdt=p_wdt,
  p_cpt=p_cpt,
  p_hpt=p_hpt
)

#### Prepare inits #############
for (c in 1:4) {
  set.seed(c)
  
  init_list <- init_f_Gaussian_no_predictor()
  
  # Define output file name
  json_file <- paste0(directory,"/sampling/inits/init_Gaussian_no_predictor_chain_", c, ".json")
  
  # Save to JSON in a Stan-readable format
  write_stan_json(init_list, json_file)
}

#### Fit model ####
mod<-cmdstan_model(stan_file=paste0(directory,"/stan_models/no_predictor/Gaussian_no_predictor.stan"),
                   dir=paste0(directory,"/sampling/c++_models"),
                   stanc_options = list("O1")
                   )

mt=12
ad=.9
iw=1000
is=1000

fit<-mod$sample(data=data,
                output_dir=paste0(directory,"/sampling/samples"),
                seed=12345,
                chains=4,
                parallel_chains = 4,
                iter_warmup =iw,
                iter_sampling= is,
                save_warmup=F,
                max_treedepth = mt,
                adapt_delta = ad,
                refresh = is/10,
                init=c(paste0(directory,"/sampling/inits/init_Gaussian_no_predictor_chain_",1:4, ".json")),
                show_messages=T,
                show_exceptions=T)

fit$save_object(paste0(directory,"/results/fits/no_predictor/Gaussian_no_predictor.rds"))

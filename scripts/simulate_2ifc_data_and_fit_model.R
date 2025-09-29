# In this script we use the posterior from the no predictor ageing and neuropathy fit
# to simulate roughly plausible 2IFC data and use this data to fit 
# a Weibull model from which we can easily derive a prior for future 2 IFC experiments
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)

#### Load fit to empirical YN data ####
fit <- readRDS("results/fits/no_predictor/Gaussian_no_predictor.rds")
draws<-fit$draws(variables = c('mu','tau'),format = 'df')

#### Simulate 2 IFC data ####
N_subj=1000
N_rep=50

stim_range=seq(.5,15,.5)
stim_range_length=length(stim_range)

y_cdt<-y_wdt<-vector(,N_subj*stim_range_length)

for(subj_idx in 1:N_subj){
  
  draw_idx<-sample(4000,1)
  
  selected_draws<-draws %>% filter(.draw==draw_idx)
  
  mu_cdt=rnorm(1,selected_draws$`mu[1]`,selected_draws$`tau[1]`)
  mu_wdt=rnorm(1,selected_draws$`mu[2]`,selected_draws$`tau[2]`)
  sigma_cdt=1/exp(rnorm(1,selected_draws$`mu[5]`,selected_draws$`tau[5]`))
  sigma_wdt=1/exp(rnorm(1,selected_draws$`mu[6]`,selected_draws$`tau[6]`))
  lambda=0.5/(1+exp(-rnorm(1,selected_draws$`mu[10]`,selected_draws$`tau[10]`)))
  
  for(stim_idx in 1:stim_range_length){
    stim=stim_range[stim_idx]
    y_cdt[stim_idx+stim_range_length*(subj_idx-1)]=rbinom(1,N_rep,0.5+(0.5-lambda)*pnorm(stim,mu_cdt,sigma_cdt))
    y_wdt[stim_idx+stim_range_length*(subj_idx-1)]=rbinom(1,N_rep,0.5+(0.5-lambda)*pnorm(stim,mu_wdt,sigma_wdt))
  }
}

n_cdt<-n_wdt<-rep(N_rep,N_subj*stim_range_length)
x_cdt<-x_wdt<-rep(stim_range,N_subj)
p_cdt<-p_wdt<-sort(rep(1:N_subj,stim_range_length))

tibble(n_cdt,y_cdt,x_cdt,n_wdt,y_wdt,x_wdt,idx=1:30000) %>% 
  pivot_longer(cols = !idx,names_to = c('var','task'),names_sep = '_') %>% 
  pivot_wider(names_from = var,values_from = value) %>% 
  group_by(task,x) %>% 
  summarise(m=median(y/n),lb=quantile(y/n,.025),ub=quantile(y/n,0.975)) %>% 
  ggplot()+
  geom_line(aes(x,y=m,color=task))+
  geom_line(aes(x,y=lb,color=task),linetype='dotted')+
  geom_line(aes(x,y=ub,color=task),linetype='dotted')+
  theme_classic()

data=list(
  P=N_subj,
  N_cdt=stim_range_length*N_subj,
  N_wdt=stim_range_length*N_subj,
  x_cdt=x_cdt,
  y_cdt=y_cdt,
  n_cdt=n_cdt,
  p_cdt=p_cdt,
  x_wdt=x_wdt,
  y_wdt=y_wdt,
  n_wdt=n_wdt,
  p_wdt=p_wdt
)

#### create init_function ####
init_f<-function(){
  init<-list(
    mu=c(
      rnorm(1,0,.5),
      rnorm(1,0,.5),
      rnorm(1,1.5,.5),
      rnorm(1,1.5,.5),
      rnorm(1,-4,1)
    ),
    tau=c(
      abs(rnorm(1,0,.5)),
      abs(rnorm(1,0,.5)),
      abs(rnorm(1,0,.5)),
      abs(rnorm(1,0,.5)),
      abs(rnorm(1,0,1))
    )
  )
  return(init)
}

#### Prepare inits #############
for (c in 1:4) {
  set.seed(c)
  
  init_list <- init_f()
  
  # Define output file name
  json_file <- paste0("sampling/inits/init_Weibull_2ifc_chain_", c, ".json")
  
  # Save to JSON in a Stan-readable format
  write_stan_json(init_list, json_file)
}

#### fit data ####
mod<-cmdstan_model(stan_file='stan_models/no_predictor/Weibull_2ifc_no_predictor.stan',
                   stanc_options = list("O1")
                   )
fit<-mod$sample(data=data,
                seed=12345,
                chains=4,
                parallel_chains = 4,
                iter_warmup =2000,
                iter_sampling= 2000,
                max_treedepth = 12,
                adapt_delta = .99,
                init=c(paste0("sampling/inits/init_Weibull_2ifc_chain_",1:4, ".json")),
                
                )

fit$save_object("results/fits/no_predictor/Weibull_2ifc_no_predictor.rds")


# Script to fit models with weakly informative priors to data acquired 
# in a previous experiment, in order to derive informative priors 
# for the ageing and neuropathy experiment
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### set-up environment
library(cmdstanr)
library(tidyverse)
options(cmdstanr_max_rows=40)

if (!dir.exists("prior_justification/fits")) {
  dir.create("prior_justification/fits", recursive = TRUE)
}

#### import and reformat data ####
d <- 
  read.csv("prior_justification/fmri_thresholding_data.csv")

d %>% 
  filter(c==1) %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE)+
  facet_wrap(p~.)
d %>% 
  filter(c==2) %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE)+
  facet_wrap(p~.)
d %>% 
  filter(c==3) %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE)+
  facet_wrap(p~.)
d %>% 
  filter(c==4) %>% 
  ggplot(aes(x=x,y=y))+
  geom_point()+
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = FALSE)+
  facet_wrap(p~.)

u<-unique(d$p)

for(idx in 1:length(u)){
  d$p[d$p==u[idx]] =idx
}

d$c<-as.factor(d$c)

x_cdt<-d$x[d$c==1]
x_wdt<-d$x[d$c==2]
x_cpt<-d$x[d$c==3]
x_hpt<-d$x[d$c==4]

y_cdt<-d$y[d$c==1]
y_wdt<-d$y[d$c==2]
y_cpt<-d$y[d$c==3]
y_hpt<-d$y[d$c==4]

p_cdt<-d$p[d$c==1]
p_wdt<-d$p[d$c==2]
p_cpt<-d$p[d$c==3]
p_hpt<-d$p[d$c==4]
  
N_cdt<-length(x_cdt)
N_wdt<-length(x_wdt)
N_cpt<-length(x_cpt)
N_hpt<-length(x_hpt)

P<-length(unique(d$p))

data<-
  list(
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
    p_hpt=p_hpt,    
    N_cdt=N_cdt,
    N_wdt=N_wdt,
    N_cpt=N_cpt,
    N_hpt=N_hpt,
    P=P
  )

#### fit and save models ####
models<-c(
  'Gaussian.stan',
  'Gaussian_with_constrained_threshold.stan',
  'Quick.stan'
)
for(m in 1:length(models)){

  m<-models[m]
  
  mod<-cmdstan_model(
    stan_file = paste('prior_justification/stan_models/',m,sep=''),
  )
  fit<-mod$sample(
    data=data,
    seed=12345,
    chains=4,
    parallel_chains=2,
    iter_warmup=2000,
    iter_sampling=2000,
    max_treedepth=12,
    adapt_delta=.95,
    refresh=100,
    init=2
  )
  fit$save_object(paste('prior_justification/fits/',str_remove(m,'.stan'),'_fit.rds',sep=''))
}
# Script to compute posterior probabilities and Bayes factors used for tests
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)
library(ggpubr)
library(loo)

rm(list=ls())
inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}
density_at_zero <- function(x) {
  dens <- density(x, from = 0, to = 0, n = 1)  # force evaluation at 0
  dens$y[1]
}

#### Load fitted models ####
fit_a<-readRDS("results/fits/Gaussian_constrained_threshold_ageing.rds")
fit_an<-readRDS("results/fits/Gaussian_constrained_threshold_ageing_and_neuropathy_NRS_NI.rds")

#### Effect of age ####
#Extract age effect coefficients from the ageing only fit
age_effect<-
  fit_a$draws(variables = paste0('mu[',11:20,']'),format = 'df') %>%
  rename(
    CD_T=`mu[11]`,
    WD_T=`mu[12]`,
    CP_T=`mu[13]`,
    HP_T=`mu[14]`,
    CD_S=`mu[15]`,
    WD_S=`mu[16]`,
    CP_S=`mu[17]`,
    HP_S=`mu[18]`,
    ALL_G=`mu[19]`,
    ALL_L=`mu[20]`
  ) %>%
  pivot_longer(
    cols = !c(.iteration,.chain,.draw),
    names_to = c("task", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(task=factor(task,c('CD','WD','CP','HP','ALL'))) %>%
  select(task,parameter,value)

#Summarize posterior
P_A<-
  age_effect %>%
  group_by(task,parameter) %>%
  summarise(
    es=sprintf("κ = %.3f [%.3f;%.3f]",mean(value),quantile(value,0.05),quantile(value,0.95)),
    p=sprintf("P(κ≤0) = %.3f",mean(value <= 0)),
    one_minus_p=sprintf("P(κ≥0) = %.3f",mean(value >= 0))
  )

#Sample from priors and compute Bayes factors
M<-8000
set.seed(12345)
age_effect_priors<-
  tibble(
    CD_T=rnorm(M,0,0.92/60),
    WD_T=rnorm(M,0,0.53/60),
    CP_T=rnorm(M,0,0.38/60),
    HP_T=rnorm(M,0,0.17/60),
    CD_S=rnorm(M,0,1.09/60),
    WD_S=rnorm(M,0,0.75/60),
    CP_S=rnorm(M,0,0.53/60),
    HP_S=rnorm(M,0,0.85/60),
    ALL_G=rnorm(M,0,0.63/60),
    ALL_L=rnorm(M,0,0.55/60)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("task", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    task=factor(task,c('CD','WD','CP','HP','ALL')),
    distribution='prior'
  )

AE<-
  age_effect %>%
  mutate(distribution='posterior') %>%
  full_join(age_effect_priors)
dens_A<-
  AE %>%
  group_by(task,parameter,distribution) %>%
  summarise(density=density_at_zero(value))
BF_A<-
  dens_A%>%
  pivot_wider(names_from = distribution,values_from = density) %>%
  mutate(BF01=sprintf("log₁₀(BF₀₁) = %.2f",log10(posterior/prior)))

#### Effect of neuropathy ####
#Extract neuropathy effect coefficients from the ageing and neuropathy fit
neuropathy_effect<-
  fit_an$draws(variables = paste0('mu[',11:20,']'),format = 'df') %>%
  rename(
    CD_T=`mu[11]`,
    WD_T=`mu[12]`,
    CP_T=`mu[13]`,
    HP_T=`mu[14]`,
    CD_S=`mu[15]`,
    WD_S=`mu[16]`,
    CP_S=`mu[17]`,
    HP_S=`mu[18]`,
    ALL_G=`mu[19]`,
    ALL_L=`mu[20]`
  ) %>%
  pivot_longer(
    cols = !c(.iteration,.chain,.draw),
    names_to = c("task", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(task=factor(task,c('CD','WD','CP','HP','ALL'))) %>%
  select(task,parameter,value)

#Summarize posterior
P_N<-
  neuropathy_effect %>%
  group_by(task,parameter) %>%
  summarise(
    es=sprintf("κ = %.3f [%.3f;%.3f]",mean(value),quantile(value,0.05),quantile(value,0.95)),
    p=sprintf("P(κ≤0) = %.3f",mean(value <= 0)),
    one_minus_p=sprintf("P(κ≥0) = %.3f",mean(value >= 0))
  )
#Sample from priors and compute Bayes factors
M<-8000
set.seed(12345)
neuropathy_effect_priors<-
  tibble(
    CD_T=rnorm(M,0,0.92),
    WD_T=rnorm(M,0,0.53),
    CP_T=rnorm(M,0,0.38),
    HP_T=rnorm(M,0,0.17),
    CD_S=rnorm(M,0,1.09),
    WD_S=rnorm(M,0,0.75),
    CP_S=rnorm(M,0,0.53),
    HP_S=rnorm(M,0,0.85),
    ALL_G=rnorm(M,0,0.63),
    ALL_L=rnorm(M,0,0.55)
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("task", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    task=factor(task,c('CD','WD','CP','HP','ALL')),
    distribution='prior'
  )

NE<-
  neuropathy_effect %>% 
  mutate(distribution='posterior') %>% 
  full_join(neuropathy_effect_priors)

dens_N<-
  NE %>%
  group_by(task,parameter,distribution) %>% 
  summarise(density=density_at_zero(value)) 

BF_N<-
  dens_N%>% 
  pivot_wider(names_from = distribution,values_from = density) %>% 
  mutate(BF01=sprintf("log₁₀(BF₀₁) = %.2f",log10(posterior/prior)))


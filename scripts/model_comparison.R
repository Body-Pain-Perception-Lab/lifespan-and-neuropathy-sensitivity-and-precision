# Script to extract LOO-CV diagnostic information and results to be reported in the paper
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(cmdstanr)
library(tidyverse)
library(loo)

rm(list=ls())
set.seed(12348)
directory=getwd()

#### Ageing ####
loo_GU<-readRDS("results/loo/ageing/MM_Gaussian_unconstrained_threshold_ageing_GLN.rds")
loo_GC<-readRDS("results/loo/ageing/MM_Gaussian_constrained_threshold_ageing_GLN.rds")
loo_QC<-readRDS("results/loo/ageing/MM_Quick_constrained_threshold_ageing_GLN.rds")

sprintf(
  "%i Pareto k's out of %i were problematic for the GU ageing model",
  sum(loo_GU$diagnostics$pareto_k>=0.7),
  length(loo_GU$diagnostics$pareto_k)
)
sprintf(
  "%i Pareto k's out of %i were problematic for the GC ageing model",
  sum(loo_GC$diagnostics$pareto_k>=0.7),
  length(loo_GC$diagnostics$pareto_k)
)
sprintf(
  "%i Pareto k's out of %i were problematic for the QC ageing model",
  sum(loo_QC$diagnostics$pareto_k>=0.7),
  length(loo_QC$diagnostics$pareto_k)
)

comparison_results_ageing<-
  loo_compare(list(GU=loo_GU,GC=loo_GC,QC=loo_QC)) %>% 
  as.data.frame() %>% 
  select(elpd_diff,se_diff) %>% 
  mutate(Z=round(elpd_diff/se_diff,2),P=round(pnorm(Z),3))
print(comparison_results_ageing)

#### Ageing and neuropathy ####
loo_GU<-readRDS("results/loo/ageing_and_neuropathy/MM_Gaussian_unconstrained_threshold_ageing_and_neuropathy_GLN_NI_H.rds")
loo_GC<-readRDS("results/loo/ageing_and_neuropathy/MM_Gaussian_constrained_threshold_ageing_and_neuropathy_GLN_NI_H.rds")
loo_QC<-readRDS("results/loo/ageing_and_neuropathy/MM_Quick_constrained_threshold_ageing_and_neuropathy_GLN_NI_H.rds")

sprintf(
  "%i Pareto k's out of %i were problematic for the GU ageing and neuropathy model",
  sum(loo_GU$diagnostics$pareto_k>=0.7),
  length(loo_GU$diagnostics$pareto_k)
)
sprintf(
  "%i Pareto k's out of %i were problematic for the GC ageing and neuropathy model",
  sum(loo_GC$diagnostics$pareto_k>=0.7),
  length(loo_GC$diagnostics$pareto_k)
)
sprintf(
  "%i Pareto k's out of %i were problematic for the QC ageing and neuropathy model",
  sum(loo_QC$diagnostics$pareto_k>=0.7),
  length(loo_QC$diagnostics$pareto_k)
)

comparison_results_ageing_and_neuropathy<-
  loo_compare(list(GU=loo_GU,GC=loo_GC,QC=loo_QC)) %>% 
  as.data.frame() %>% 
  select(elpd_diff,se_diff) %>% 
  mutate(Z=round(elpd_diff/se_diff,2),P=round(pnorm(Z),3))
print(comparison_results_ageing_and_neuropathy)

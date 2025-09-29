# Script to test model recovery and generate model recovery fiGUre
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up the environement ####
library(cmdstanr)
library(tidyverse)
library(loo)

#### Extract and compare LOO-CV objects
model_recovery<-list()
model_recovery_df<-
  tibble(
    simulation=sort(rep(1:12,3)),
    gen=rep(c('GU','GC','QC'),12),
    win=vector(,12*3)
  )

for (idx in c(1:12)){
  loo_GU_GU<-readRDS(paste('parameter_recovery/loo/GU_GU_',idx,'.rds',sep=''))
  if(sum(loo_GU_GU$diagnostics$pareto_k>0.7)>length(loo_GU_GU$diagnostics$pareto_k)/50){
    loo_GU_GU$elpd_loo=-10^10
  }
  
  loo_GU_GC<-readRDS(paste('parameter_recovery/loo/GU_GC_',idx,'.rds',sep=''))
  if(sum(loo_GU_GC$diagnostics$pareto_k>0.7)>length(loo_GU_GC$diagnostics$pareto_k)/50){
    loo_GU_GC$elpd_loo=-10^10
  }
  
  loo_GU_QC<-readRDS(paste('parameter_recovery/loo/GU_QC_',idx,'.rds',sep=''))
  if(sum(loo_GU_QC$diagnostics$pareto_k>0.7)>length(loo_GU_QC$diagnostics$pareto_k)/50){
    loo_GU_QC$elpd_loo=-10^10
  }
  
  
  loo_GC_GU<-readRDS(paste('parameter_recovery/loo/GC_GU_',idx,'.rds',sep=''))
  if(sum(loo_GC_GU$diagnostics$pareto_k>0.7)>length(loo_GC_GU$diagnostics$pareto_k)/50){
    loo_GC_GU$elpd_loo=-10^10
  }
  
  loo_GC_GC<-readRDS(paste('parameter_recovery/loo/GC_GC_',idx,'.rds',sep=''))
  if(sum(loo_GC_GC$diagnostics$pareto_k>0.7)>length(loo_GC_GC$diagnostics$pareto_k)/50){
    loo_GC_GC$elpd_loo=-10^10
  }
  
  loo_GC_QC<-readRDS(paste('parameter_recovery/loo/GC_QC_',idx,'.rds',sep=''))
  if(sum(loo_GC_QC$diagnostics$pareto_k>0.7)>length(loo_GC_QC$diagnostics$pareto_k)/50){
    loo_GC_QC$elpd_loo=-10^10
  }
  
  
  loo_QC_GU<-readRDS(paste('parameter_recovery/loo/QC_GU_',idx,'.rds',sep=''))
  if(sum(loo_QC_GU$diagnostics$pareto_k>0.7)>length(loo_QC_GU$diagnostics$pareto_k)/50){
    loo_QC_GU$elpd_loo=-10^10
  }
  
  loo_QC_GC<-readRDS(paste('parameter_recovery/loo/QC_GC_',idx,'.rds',sep=''))
  if(sum(loo_QC_GC$diagnostics$pareto_k>0.7)>length(loo_QC_GC$diagnostics$pareto_k)/50){
    loo_QC_GC$elpd_loo=-10^10
  }
  
  loo_QC_QC<-readRDS(paste('parameter_recovery/loo/QC_QC_',idx,'.rds',sep=''))
  if(sum(loo_QC_QC$diagnostics$pareto_k>0.7)>length(loo_QC_QC$diagnostics$pareto_k)/50){
    loo_QC_QC$elpd_loo=-10^10
  }

  model_comp<-list(
    GU=loo_compare(list(GU=loo_GU_GU,GC=loo_GU_GC,QC=loo_GU_QC)),
    GC=loo_compare(list(GU=loo_GC_GU,GC=loo_GC_GC,QC=loo_GC_QC)),
    QC=loo_compare(list(GU=loo_QC_GU,GC=loo_QC_GC,QC=loo_QC_QC))
  )
  
  model_recovery[[idx]]=list(
    model_comp=model_comp
  )
  model_recovery_df$win[(model_recovery_df$gen=='GU')&(model_recovery_df$simulation==idx)]=rownames(model_comp[[1]])[1]
  model_recovery_df$win[(model_recovery_df$gen=='GC')&(model_recovery_df$simulation==idx)]=rownames(model_comp[[2]])[1]
  model_recovery_df$win[(model_recovery_df$gen=='QC')&(model_recovery_df$simulation==idx)]=rownames(model_comp[[3]])[1]
}

model_recovery_df<-
  model_recovery_df %>% 
  group_by(win,gen) %>% 
  summarise(n=sum(!is.na(simulation)))


model_recovery_df$gen[model_recovery_df$gen=='GU']='1b'
model_recovery_df$gen[model_recovery_df$gen=='GC']='2b'
model_recovery_df$gen[model_recovery_df$gen=='QC']='3b'
model_recovery_df$win[model_recovery_df$win=='GU']='1b'
model_recovery_df$win[model_recovery_df$win=='GC']='2b'
model_recovery_df$win[model_recovery_df$win=='QC']='3b'

model_recovery_df %>%
  ggplot()+
  geom_tile(aes(x=gen,y=win,fill=n/15))+
  geom_text(aes(x=gen,y=win,label=paste0(round(n/15,2),'\n [',round(qbeta(.025,n+1,(15-n)+1),2),' - ',round(qbeta(.975,n+1,(15-n)+1),2),']')),color='grey')+
  xlab('Generative model')+
  ylab('Winning model')+
  theme_minimal()+
  theme(legend.position = 'none')






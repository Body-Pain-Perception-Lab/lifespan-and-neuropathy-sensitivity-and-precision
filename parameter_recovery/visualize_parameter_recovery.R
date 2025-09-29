# Script to test parameter recovery and generate model recovery figure
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set up environment ####
library(cmdstanr)
library(tidyverse)
library(ggpubr)

########## Qc ##################
all_datasets <- read.csv("Guc_dataset_info.csv")
all_datasets<-data.matrix(all_datasets,rownames.force = F)

recovery<-tibble(mu_mean=rep(NaN,40*15))
recovery$mu_lb<-recovery$mu_mean
recovery$mu_ub<-recovery$mu_mean
recovery$mu<-recovery$mu_mean
recovery$tau<-recovery$mu_mean
recovery$tau_mean<-recovery$mu_mean
recovery$tau_lb<-recovery$mu_mean
recovery$tau_ub<-recovery$mu_mean

for (idx in c(1:15)){
  fit<-readRDS(paste('fits/GU_GU_',idx,'.rds',sep=''))
  
  s<-fit$summary(variables = c('mu','tau'))
  
  mu_mean<-s$mean[1:40]
  tau_mean<-s$mean[41:60]
  mu_lb<-s$q5[1:40]
  tau_lb<-s$q5[41:60]
  mu_ub<-s$q95[1:40]
  tau_ub<-s$q95[41:60]
  
  mu<-all_datasets[idx,2:41]
  tau<-all_datasets[idx,42:61]
  
  recovery$mu_mean[(1:40)+(idx-1)*40]<-mu_mean
  recovery$mu_lb[(1:40)+(idx-1)*40]<-mu_lb
  recovery$mu_ub[(1:40)+(idx-1)*40]<-mu_ub
  recovery$mu[(1:40)+(idx-1)*40]<-mu
  
  recovery$tau_mean[(1:20)+(idx-1)*40]<-tau_mean
  recovery$tau_lb[(1:20)+(idx-1)*40]<-tau_lb
  recovery$tau_ub[(1:20)+(idx-1)*40]<-tau_ub 
  recovery$tau[(1:20)+(idx-1)*40]<-tau
  
}
recovery_GU<-recovery

recovery_GU %>%
  ggplot()+
  geom_linerange(aes(x=mu,y=mu_mean,ymin=mu_lb,ymax=mu_ub),alpha=.5)+
  geom_abline(linetype='dotted')+
  ylab('mu estimate')+
  xlab('mu')+
  theme_classic()


recovery_GU %>%
  ggplot()+
  geom_point(aes(x=log(abs(mu)),y=(mu_mean-mu)/(mu_ub-mu_lb)),alpha=.5,size=.5)+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=-.5),linetype='dashed')+
  geom_hline(aes(yintercept=.5),linetype='dashed')+
  ylab('standardized deviation from ground truth')+
  xlab('log(abs(mu))')+
  theme_classic()

recovery_GU %>% 
  ggplot()+
  geom_linerange(aes(x=log(tau),y=log(tau_mean),ymin=log(tau_lb),ymax=log(tau_ub)),alpha=.5)+
  geom_abline(linetype='dotted')+
  ylab('log(tau) estimate')+
  xlab('log(tau)')+
  theme_classic()

recovery_GU %>%
  ggplot()+
  geom_point(aes(x=log(tau),y=(tau_mean-tau)/(tau_ub-tau_lb)),alpha=.5,size=.5)+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=-.5),linetype='dashed')+
  geom_hline(aes(yintercept=.5),linetype='dashed')+
  ylab('standardized deviation from ground truth')+
  xlab('log(tau)')+
  theme_classic()
########## Gc ##################
all_datasets <- read.csv("Gc_dataset_info.csv")
all_datasets<-data.matrix(all_datasets,rownames.force = F)

recovery<-tibble(mu_mean=rep(NaN,40*15))
recovery$mu_lb<-recovery$mu_mean
recovery$mu_ub<-recovery$mu_mean
recovery$mu<-recovery$mu_mean
recovery$tau<-recovery$mu_mean
recovery$tau_mean<-recovery$mu_mean
recovery$tau_lb<-recovery$mu_mean
recovery$tau_ub<-recovery$mu_mean

for (idx in c(1:15)){
  fit<-readRDS(paste('fits/GC_GC_',idx,'.rds',sep=''))

  s<-fit$summary(variables = c('mu','tau'))
  
  mu_mean<-s$mean[1:40]
  tau_mean<-s$mean[41:60]
  mu_lb<-s$q5[1:40]
  tau_lb<-s$q5[41:60]
  mu_ub<-s$q95[1:40]
  tau_ub<-s$q95[41:60]
  
  mu<-all_datasets[idx,2:41]
  tau<-all_datasets[idx,42:61]
  
  recovery$mu_mean[(1:40)+(idx-1)*40]<-mu_mean
  recovery$mu_lb[(1:40)+(idx-1)*40]<-mu_lb
  recovery$mu_ub[(1:40)+(idx-1)*40]<-mu_ub
  recovery$mu[(1:40)+(idx-1)*40]<-mu
  
  recovery$tau_mean[(1:20)+(idx-1)*40]<-tau_mean
  recovery$tau_lb[(1:20)+(idx-1)*40]<-tau_lb
  recovery$tau_ub[(1:20)+(idx-1)*40]<-tau_ub 
  recovery$tau[(1:20)+(idx-1)*40]<-tau
  
}
recovery_GC<-recovery


recovery_GC %>%
  ggplot()+
  geom_linerange(aes(x=mu,y=mu_mean,ymin=mu_lb,ymax=mu_ub),alpha=.5)+
  geom_abline(linetype='dotted')+
  ylab('mu estimate')+
  xlab('mu')+
  theme_classic()


recovery_GC %>%
  ggplot()+
  geom_point(aes(x=log(abs(mu)),y=(mu_mean-mu)/(mu_ub-mu_lb)),alpha=.5,size=.5)+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=-.5),linetype='dashed')+
  geom_hline(aes(yintercept=.5),linetype='dashed')+
  ylab('standardized deviation from ground truth')+
  xlab('log(abs(mu))')+
  theme_classic()

recovery_GC %>% 
  ggplot()+
  geom_linerange(aes(x=log(tau),y=log(tau_mean),ymin=log(tau_lb),ymax=log(tau_ub)),alpha=.5)+
  geom_abline(linetype='dotted')+
  ylab('log(tau) estimate')+
  xlab('log(tau)')+
  theme_classic()

recovery_GC %>%
  ggplot()+
  geom_point(aes(x=log(tau),y=(tau_mean-tau)/(tau_ub-tau_lb)),alpha=.5,size=.5)+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=-.5),linetype='dashed')+
  geom_hline(aes(yintercept=.5),linetype='dashed')+
  ylab('standardized deviation from ground truth')+
  xlab('log(tau)')+
  theme_classic()
########## Qc ##################
all_datasets <- read.csv("Qc_dataset_info.csv")
all_datasets<-data.matrix(all_datasets,rownames.force = F)

recovery<-tibble(mu_mean=rep(NaN,40*15))
recovery$mu_lb<-recovery$mu_mean
recovery$mu_ub<-recovery$mu_mean
recovery$mu<-recovery$mu_mean
recovery$tau<-recovery$mu_mean
recovery$tau_mean<-recovery$mu_mean
recovery$tau_lb<-recovery$mu_mean
recovery$tau_ub<-recovery$mu_mean

for (idx in c(1:15)){
  fit<-readRDS(paste('fits/QC_QC_',idx,'.rds',sep=''))

  s<-fit$summary(variables = c('mu','tau'))
  
  mu_mean<-s$mean[1:40]
  tau_mean<-s$mean[41:60]
  mu_lb<-s$q5[1:40]
  tau_lb<-s$q5[41:60]
  mu_ub<-s$q95[1:40]
  tau_ub<-s$q95[41:60]
  
  mu<-all_datasets[idx,2:41]
  tau<-all_datasets[idx,42:61]
  
  recovery$mu_mean[(1:40)+(idx-1)*40]<-mu_mean
  recovery$mu_lb[(1:40)+(idx-1)*40]<-mu_lb
  recovery$mu_ub[(1:40)+(idx-1)*40]<-mu_ub
  recovery$mu[(1:40)+(idx-1)*40]<-mu
  
  recovery$tau_mean[(1:20)+(idx-1)*40]<-tau_mean
  recovery$tau_lb[(1:20)+(idx-1)*40]<-tau_lb
  recovery$tau_ub[(1:20)+(idx-1)*40]<-tau_ub 
  recovery$tau[(1:20)+(idx-1)*40]<-tau
  
}
recovery_QC<-recovery

recovery_QC %>%
  ggplot()+
  geom_linerange(aes(x=mu,y=mu_mean,ymin=mu_lb,ymax=mu_ub),alpha=.5)+
  geom_abline(linetype='dotted')+
  ylab('mu estimate')+
  xlab('mu')+
  theme_classic()


recovery_QC %>%
  ggplot()+
  geom_point(aes(x=log(abs(mu)),y=(mu_mean-mu)/(mu_ub-mu_lb)),alpha=.5,size=.5)+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=-.5),linetype='dashed')+
  geom_hline(aes(yintercept=.5),linetype='dashed')+
  ylab('standardized deviation from ground truth')+
  xlab('log(abs(mu))')+
  theme_classic()

recovery_QC %>% 
  ggplot()+
  geom_linerange(aes(x=log(tau),y=log(tau_mean),ymin=log(tau_lb),ymax=log(tau_ub)),alpha=.5)+
  geom_abline(linetype='dotted')+
  ylab('log(tau) estimate')+
  xlab('log(tau)')+
  theme_classic()

recovery_QC %>%
  ggplot()+
  geom_point(aes(x=log(tau),y=(tau_mean-tau)/(tau_ub-tau_lb)),alpha=.5,size=.5)+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=-.5),linetype='dashed')+
  geom_hline(aes(yintercept=.5),linetype='dashed')+
  ylab('standardized deviation from ground truth')+
  xlab('log(tau)')+
  theme_classic()

############# Create plot ############
mu_GU<-
  recovery_GU %>%
  ggplot()+
  geom_abline(color='gray')+
  geom_linerange(aes(x=mu,y=mu_mean,ymin=mu_lb,ymax=mu_ub),alpha=.2)+
  ylab('estimate')+
  xlab('mu')+
  theme_classic()

mu_GC<-
  recovery_GC %>%
  ggplot()+
  geom_abline(color='gray')+
  geom_linerange(aes(x=mu,y=mu_mean,ymin=mu_lb,ymax=mu_ub),alpha=.2)+
  ylab('estimate')+
  xlab('mu')+
  theme_classic()

mu_QC<-
  recovery_QC %>%
  ggplot()+
  geom_abline(color='gray')+
  geom_linerange(aes(x=mu,y=mu_mean,ymin=mu_lb,ymax=mu_ub),alpha=.2)+
  ylab('estimate')+
  xlab('mu')+
  theme_classic()

tau_GC<-
  recovery_GC %>% 
  ggplot()+
  geom_abline(color='gray')+
  geom_linerange(aes(x=tau,y=(tau_mean),ymin=(tau_lb),ymax=(tau_ub)),alpha=.2)+
  ylab('')+
  xlab('tau')+
  theme_classic()

tau_GU<-
  recovery_GU %>% 
  ggplot()+
  geom_abline(color='gray')+
  geom_linerange(aes(x=tau,y=(tau_mean),ymin=(tau_lb),ymax=(tau_ub)),alpha=.2)+
  ylab('')+
  xlab('tau')+
  theme_classic()

tau_QC<-
  recovery_QC %>% 
  ggplot()+
  geom_abline(color='gray')+
  geom_linerange(aes(x=tau,y=(tau_mean),ymin=(tau_lb),ymax=(tau_ub)),alpha=.2)+
  ylab('')+
  xlab('tau')+
  theme_classic()

ggarrange(mu_GU,tau_GU,mu_GC,tau_GC,mu_QC,tau_QC,ncol = 2,nrow = 3)




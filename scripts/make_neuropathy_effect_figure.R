# Script to generate figure 5
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)
library(ggpubr)

rm(list=ls())
inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}

#### Load fitted models ####
fit_an<-readRDS("results/fits/Gaussian_constrained_threshold_ageing_and_neuropathy_NRS_NI.rds")

#### Effect of neuropathy ####
#Extract participant ages
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

age<-tibble(
  age=participant_info$age,
  status=participant_info$ID>5000,
)
age$participant<-1:length(age$status)

#Extract participant-level parameter estimates
ppt_params_N<-
  fit_an$draws(c('alpha','beta'),format='df') %>% 
  pivot_longer(!c(.iteration,.draw,.chain)) %>% 
  group_by(name) %>% 
  summarise(
    m=mean(value),
    lb=quantile(value,.025),
    ub=quantile(value,.975)
  ) %>%  
  extract(
    name,
    into = c("parameter", "condition", "participant"),
    regex = "([A-Za-z]+)\\[(\\d+),\\s*(\\d+)\\]"
  ) %>% 
  mutate(
    condition=factor(as.numeric(condition),1:4,c('CD','WD','CP','HP')),
    participant=as.numeric(participant),
  ) %>% 
  full_join(age) %>% 
  mutate(status=factor(as.numeric(status),0:1,c('Control','Patient')))

#Estimate and summarize group mean parameters as a function of age and neuropathy
predictors<-expand_grid(.draw=1:8000,age=seq(15,85,.5),status=c(F,T))

parameters_an<-fit_an$draws(variables = 'mu',format='df')%>%
  rename(
    cd_ti=`mu[1]`,
    wd_ti=`mu[2]`,
    cp_ti=`mu[3]`,
    hp_ti=`mu[4]`,
    cd_si=`mu[5]`,
    wd_si=`mu[6]`,
    cp_si=`mu[7]`,
    hp_si=`mu[8]`,
    cd_tsn=`mu[11]`,
    wd_tsn=`mu[12]`,
    cp_tsn=`mu[13]`,
    hp_tsn=`mu[14]`,
    cd_ssn=`mu[15]`,
    wd_ssn=`mu[16]`,
    cp_ssn=`mu[17]`,
    hp_ssn=`mu[18]`,
    cd_tsa=`mu[21]`,
    wd_tsa=`mu[22]`,
    cp_tsa=`mu[23]`,
    hp_tsa=`mu[24]`,
    cd_ssa=`mu[25]`,
    wd_ssa=`mu[26]`,
    cp_ssa=`mu[27]`,
    hp_ssa=`mu[28]`,
  ) %>%
  full_join(predictors) %>%
  mutate(
    ra=age-20,
    alpha_1=exp(cd_ti+cd_tsa*ra+cd_tsn*status),
    alpha_2=exp(wd_ti+wd_tsa*ra+wd_tsn*status),
    alpha_3=exp(cp_ti+cp_tsa*ra+cp_tsn*status),
    alpha_4=exp(hp_ti+hp_tsa*ra+hp_tsn*status),
    beta_1=exp(cd_si+cd_ssa*ra+cd_ssn*status),
    beta_2=exp(wd_si+wd_ssa*ra+wd_ssn*status),
    beta_3=exp(cp_si+cp_ssa*ra+cp_ssn*status),
    beta_4=exp(hp_si+hp_ssa*ra+hp_ssn*status)
  ) %>%
  pivot_longer(
    cols = c(alpha_1,alpha_2,alpha_3,alpha_4,beta_1,beta_2,beta_3,beta_4)
  ) %>%
  select(age,status,name,value) %>% 
  separate_wider_delim(
    cols=name,
    delim="_",
    names = c("parameter", "condition")
  ) %>% 
  mutate(
    condition=factor(as.numeric(condition),1:4,c('CD','WD','CP','HP')),
    status=factor(as.numeric(status),0:1,c('Control','Patient'))
  ) %>% 
  group_by(age,status,parameter,condition) %>% 
  summarise(
    m=mean(value),
    lb=quantile(value,.025),
    ub=quantile(value,.975)
  )

#Estimate and summarize PFs
predictors<-expand_grid(x=seq(0.1,30,.1),.draw=1:8000,age=50,status=c(0,1))

pf_neuropathy<-fit_an$draws(variables = 'mu',format='df')%>%
  rename(
    cd_ti=`mu[1]`,
    wd_ti=`mu[2]`,
    cp_ti=`mu[3]`,
    hp_ti=`mu[4]`,
    cd_si=`mu[5]`,
    wd_si=`mu[6]`,
    cp_si=`mu[7]`,
    hp_si=`mu[8]`,
    g_i=`mu[9]`,
    l_i=`mu[10]`,
    cd_ts=`mu[11]`,
    wd_ts=`mu[12]`,
    cp_ts=`mu[13]`,
    hp_ts=`mu[14]`,
    cd_ss=`mu[15]`,
    wd_ss=`mu[16]`,
    cp_ss=`mu[17]`,
    hp_ss=`mu[18]`,
    g_s=`mu[19]`,
    l_s=`mu[20]`,
    cd_ta=`mu[21]`,
    wd_ta=`mu[22]`,
    cp_ta=`mu[23]`,
    hp_ta=`mu[24]`,
    cd_sa=`mu[25]`,
    wd_sa=`mu[26]`,
    cp_sa=`mu[27]`,
    hp_sa=`mu[28]`,
    g_a=`mu[29]`,
    l_a=`mu[30]`
  ) %>% 
  full_join(predictors) %>% 
  mutate(
    ra=age-20,
    guess=.5*inv_logit(g_i+g_a*ra+g_s*status),
    lapse=.5*inv_logit(l_i+l_a*ra+l_s*status),
    cdt=exp(cd_ti+cd_ta*ra+cd_ts*status),
    wdt=exp(wd_ti+wd_ta*ra+wd_ts*status),
    cpt=exp(cp_ti+cp_ta*ra+cp_ts*status),
    hpt=exp(hp_ti+hp_ta*ra+hp_ts*status),
    cds=exp(cd_si+cd_sa*ra+cd_ss*status),
    wds=exp(wd_si+wd_sa*ra+wd_ss*status),
    cps=exp(cp_si+cp_sa*ra+cp_ss*status),
    hps=exp(hp_si+hp_sa*ra+hp_ss*status),
    CD_p=guess+(1-guess-lapse)*pnorm(x,cdt,1/cds),
    WD_p=guess+(1-guess-lapse)*pnorm(x,wdt,1/wds),
    CP_p=guess+(1-guess-lapse)*pnorm(x,cpt,1/cps),
    HP_p=guess+(1-guess-lapse)*pnorm(x,hpt,1/hps),
  ) %>%
  pivot_longer(
    cols = ends_with("_p"),
    names_to = "task",
    names_pattern = "(.*)_p",
    values_to = "p"
  ) %>% 
  mutate(task=factor(task,c('CD','WD','CP','HP'))) %>% 
  filter((task=='CD'&x<15)|(task=='WD'&x<15)|task=='CP'|(task=='HP'&x<20)) %>% 
  group_by(x,status,task) %>% 
  summarise(
    m=median(p),
    lb=quantile(p,0.025),
    ub=quantile(p,0.975)
  ) %>% 
  mutate(status_st=if_else(status==0, "Control", "Patient"))

pf_neuropathy$x[pf_neuropathy$task=='CD']<-30-pf_neuropathy$x[pf_neuropathy$task=='CD']
pf_neuropathy$x[pf_neuropathy$task=='WD']<-30+pf_neuropathy$x[pf_neuropathy$task=='WD']
pf_neuropathy$x[pf_neuropathy$task=='CP']<-30-pf_neuropathy$x[pf_neuropathy$task=='CP']
pf_neuropathy$x[pf_neuropathy$task=='HP']<-30+pf_neuropathy$x[pf_neuropathy$task=='HP']

#Make effect of neuropathy on threshold plot
NTP <-
  ggplot() +
  geom_line(
    data=parameters_an %>%filter(parameter == 'alpha'),
    aes(x=age,color=condition,y=m,alpha=status)
  )+
  geom_ribbon(
    data=parameters_an %>%filter(parameter == 'alpha'),
    aes(x=age,ymin=lb,ymax=ub,alpha=status,color=condition),
  )+
  geom_pointrange(
    data=ppt_params_N %>%filter(parameter == 'alpha'), 
    aes(x=age,color=condition,y=m,ymin=lb,ymax=ub,alpha=status),
    size=.1,
    linewidth=.1
  )+
  scale_color_manual(
    values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('CD', 'WD', 'CP', 'HP')
  ) +  
  guides(color='none',alpha='none',linetype='none',shape='none')+
  labs(
    title='Effects of neuropathy on the thresholds',
    x = 'Age (year)',
    y = 'Threshold (°C)',
    color = 'Modality'
  ) +
  scale_alpha_manual(values=c(.2,.7))+
  facet_wrap(~condition, ncol = 4,scales = 'free') +
  theme_classic() +
  theme(strip.text.x = element_blank())+
  scale_y_log10()

#Make alternative effect of neuropathy on slope plot
NSP <-
  ggplot() +
  geom_line(
    data=parameters_an %>%filter(parameter == 'beta'),
    aes(x=age,color=condition,y=m,alpha=status)
  )+
  geom_ribbon(
    data=parameters_an %>%filter(parameter == 'beta'),
    aes(x=age,ymin=lb,ymax=ub,color=condition,alpha=status),
  )+
  geom_pointrange(
    data=ppt_params_N %>%filter(parameter == 'beta'), 
    aes(x=age,color=condition,y=m,ymin=lb,ymax=ub,alpha=status),
    size=.1,
    linewidth=.1
  )+
  scale_color_manual(
    values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('CD', 'WD', 'CP', 'HP')
  ) +  
  guides(color='none',linetype='none',shape='none')+
  labs(
    title='Effects of neuropathy on the slopes',
    x = 'Age (year)',
    y = 'Slope (°C⁻¹)',
    alpha = 'Status'
  ) +
  scale_alpha_manual(values=c(.2,.7))+
  facet_wrap(~condition, ncol = 4,scales = 'free') +
  theme_classic() +
  theme(
    strip.text.x = element_blank(),
    legend.position = 'bottom')+
  scale_y_log10()

#Make group mean PF plot
PFNP<-
  pf_neuropathy %>%
  ggplot()+
  geom_vline(xintercept = 30, linetype = 'dashed', color = 'grey') +
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub,fill=task,alpha=status_st,group=interaction(task,status)),alpha=.5)+
  geom_line(aes(x=x,y=m,color=task))+
  theme_classic()+
    ggh4x::facet_grid2(
      rows = vars(status_st),
      cols = vars(task),
      scales = "free",
      axes = "x"  # show x-axis line everywhere, breaks only on bottom
    ) +   xlab('Temperature (°C)')+
  ylab('P(response="yes")')+
  scale_color_manual(
      values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
      labels = c('CD', 'WD', 'CP', 'HP')
    ) +  
  scale_fill_manual(
      values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
      labels = c('CD', 'WD', 'CP', 'HP')
    ) +  
  labs(
    fill='Task',
    color='Task',
    title='Group mean psychometric functions at age 50 for patients and controls'
  )+
  theme(  
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_line(),
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    legend.position = 'bottom'
  )+
  scale_y_continuous(breaks = c(0,.5,1))

#Combine and save
combinedN <- ggarrange(NTP, NSP, PFNP, ncol = 1,heights = c(1,1.2,2))

ggsave(
  'plots/neuropathy_effect_plot.jpeg',
  plot = combinedN,
  width = 18,
  height = 24,
  dpi = 300,
  units = 'cm'
  )

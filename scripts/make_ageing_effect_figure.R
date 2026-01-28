# Script to generate figure 4
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(ggpubr)

rm(list=ls())
inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}

#### Load fitted model ####
fit_a<-readRDS("results/fits/Gaussian_constrained_threshold_ageing.rds")

#### Effect of age ####
#Extract participant ages
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

age<-tibble(
  age=participant_info$age,
  participant=1:75
)

#Extract participant-level parameter estimates
ppt_params<-
  fit_a$draws(c('alpha','beta'),format='df') %>% 
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
    participant=as.numeric(participant)
  ) %>% 
  full_join(age)

#Estimate and summarize group mean parameters as a function of age
predictors<-expand_grid(.draw=1:8000,age=seq(15,85,.5))

parameters_ageing<-fit_a$draws(variables = 'mu',format='df')%>%
  rename(
    cd_ti=`mu[1]`,
    wd_ti=`mu[2]`,
    cp_ti=`mu[3]`,
    hp_ti=`mu[4]`,
    cd_si=`mu[5]`,
    wd_si=`mu[6]`,
    cp_si=`mu[7]`,
    hp_si=`mu[8]`,
    cd_ts=`mu[11]`,
    wd_ts=`mu[12]`,
    cp_ts=`mu[13]`,
    hp_ts=`mu[14]`,
    cd_ss=`mu[15]`,
    wd_ss=`mu[16]`,
    cp_ss=`mu[17]`,
    hp_ss=`mu[18]`,
  ) %>%
  full_join(predictors) %>%
  mutate(
    ra=age-20,
    alpha_1=exp(cd_ti+cd_ts*ra),
    alpha_2=exp(wd_ti+wd_ts*ra),
    alpha_3=exp(cp_ti+cp_ts*ra),
    alpha_4=exp(hp_ti+hp_ts*ra),
    beta_1=exp(cd_si+cd_ss*ra),
    beta_2=exp(wd_si+wd_ss*ra),
    beta_3=exp(cp_si+cp_ss*ra),
    beta_4=exp(hp_si+hp_ss*ra)
  ) %>%
  pivot_longer(
    cols = c(alpha_1,alpha_2,alpha_3,alpha_4,beta_1,beta_2,beta_3,beta_4)
  ) %>%
  select(age,name,value) %>% 
  separate_wider_delim(
    cols=name,
    delim="_",
    names = c("parameter", "condition")
  ) %>% 
  mutate(
    condition=factor(as.numeric(condition),1:4,c('CD','WD','CP','HP')),
  ) %>% 
  group_by(age,parameter,condition) %>% 
  summarise(
    m=mean(value),
    lb=quantile(value,.025),
    ub=quantile(value,.975)
  )

#Estimate and summarize group mean PFs
predictors<-expand_grid(x=seq(0.1,30,.1),.draw=1:8000,age=c(20,80))

pf_ageing<-fit_a$draws(variables = 'mu',format='df')%>%
  rename(
    cd_ti=`mu[1]`,
    wd_ti=`mu[2]`,
    cp_ti=`mu[3]`,
    hp_ti=`mu[4]`,
    cd_si=`mu[5]`,
    wd_si=`mu[6]`,
    cp_si=`mu[7]`,
    hp_si=`mu[8]`,
    g_si=`mu[9]`,
    l_si=`mu[10]`,
    cd_ts=`mu[11]`,
    wd_ts=`mu[12]`,
    cp_ts=`mu[13]`,
    hp_ts=`mu[14]`,
    cd_ss=`mu[15]`,
    wd_ss=`mu[16]`,
    cp_ss=`mu[17]`,
    hp_ss=`mu[18]`,
    g_ss=`mu[19]`,
    l_ss=`mu[20]`,
  ) %>%
  full_join(predictors) %>%
  mutate(
    ra=age-20,
    guess=0.5*inv_logit(g_si+g_ss*ra),
    lapse=0.5*inv_logit(l_si+l_ss*ra),
    cdt=exp(cd_ti+cd_ts*ra),
    wdt=exp(wd_ti+wd_ts*ra),
    cpt=exp(cp_ti+cp_ts*ra),
    hpt=exp(hp_ti+hp_ts*ra),
    cds=exp(cd_si+cd_ss*ra),
    wds=exp(wd_si+wd_ss*ra),
    cps=exp(cp_si+cp_ss*ra),
    hps=exp(hp_si+hp_ss*ra),
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
  filter((task=='CD'&x<15)|(task=='WD'&x<15)|task=='CP'|(task=='HP'&x<20))%>%
  group_by(x,age,task) %>%
  summarise(
    m=median(p),
    lb=quantile(p,0.025),
    ub=quantile(p,0.975)
  )

pf_ageing$x[pf_ageing$task=='CD']<-30-pf_ageing$x[pf_ageing$task=='CD']
pf_ageing$x[pf_ageing$task=='WD']<-30+pf_ageing$x[pf_ageing$task=='WD']
pf_ageing$x[pf_ageing$task=='CP']<-30-pf_ageing$x[pf_ageing$task=='CP']
pf_ageing$x[pf_ageing$task=='HP']<-30+pf_ageing$x[pf_ageing$task=='HP']

#Make effect of age on threshold plot
ATP <-
  ggplot() +
  geom_line(
      data=parameters_ageing %>%filter(parameter == 'alpha'), 
      aes(x=age,color=condition,y=m)
    )+
  geom_ribbon(
      data=parameters_ageing %>%filter(parameter == 'alpha'), 
      aes(x=age,color=condition,ymin=lb,ymax=ub),
      alpha=.1
    )+
  geom_pointrange(
    data=ppt_params %>%filter(parameter == 'alpha'), 
    aes(x=age,color=condition,y=m,ymin=lb,ymax=ub),
    alpha=.5,
    size=.1,
    linewidth=.1
  )+
  scale_color_manual(
    values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('CD', 'WD', 'CP', 'HP')
  ) +
  guides(color='none',linetype='none',shape='none')+
  labs(
    title='Effect of age on the thresholds of healthy volunteers',
    x = 'Age (year)',
    y = 'Threshold (°C)',
    color = 'Modality'
  ) +
  facet_wrap(~condition, ncol = 4,scales = 'free') +
  theme_classic() +
  theme(strip.text.x = element_blank())+
  scale_y_log10()

#Make alternative effect of age on slope plot
ASP <-
  ggplot() +
  geom_line(
    data=parameters_ageing %>%filter(parameter == 'beta'),
    aes(x=age,color=condition,y=m)
  )+
  geom_ribbon(
    data=parameters_ageing %>%filter(parameter == 'beta'),
    aes(x=age,color=condition,ymin=lb,ymax=ub),
    alpha=.1
  )+
  geom_pointrange(
    data=ppt_params %>%filter(parameter == 'beta'), 
    aes(x=age,color=condition,y=m,ymin=lb,ymax=ub),
    alpha=.5,
    size=.1,
    linewidth=.1
  )+
  scale_color_manual(
    values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('CD', 'WD', 'CP', 'HP')
  ) +
  labs(
    title='Effect of age on the slopes of healthy volunteers',
    x = 'Age (year)',
    y = 'Slope (°C⁻¹)',
    color = 'Modality'
  ) +
  guides(color='none')+
  facet_wrap(~condition, ncol = 4,scales = 'free') +
  theme_classic() +
  theme(
    strip.text.x = element_blank()
    )+
  scale_y_log10()

#Make group mean PF plot
PFAP<-pf_ageing %>%
  ggplot()+
  geom_vline(xintercept = 30, linetype = 'dashed', color = 'grey') +
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub,fill=task),alpha=0.5)+
  geom_line(aes(x=x,y=m,color=task))+
  theme_classic()+
  ggh4x::facet_grid2(
    rows = vars(age),
    cols = vars(task),
    scales = "free",
    axes = "x"  # show x-axis line everywhere, breaks only on bottom
  ) +   
  scale_color_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  scale_fill_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  xlab('Temperature (°C)')+
  ylab('P(response="yes")')+
  labs(
    fill='Modality',
    color='Modality',
    title='Group mean psychometric functions of healthy volunteers at age 20 and 80'
  )+
  theme(
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    legend.position = 'bottom'
  )+
  scale_y_continuous(breaks = c(0,.5,1))

#Combine and save
combinedA <- ggarrange(ATP, ASP, PFAP, ncol = 1,heights = c(1,1,2))

ggsave(
  'plots/ageing_effect_plot.jpeg',
  plot = combinedA,
  width = 18,
  height = 24,
  dpi = 300,
  units = 'cm'
  )

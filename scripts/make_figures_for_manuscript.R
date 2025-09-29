# Script to generate all the figures in the manuscript's main text 
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
set.seed(12345)

#### Model space illustration ####
#Estimate and plot the two types of PF
PFs<-
  tibble(x=seq(0.01,2,.01)) %>% 
  mutate(
    Gaussian=.02+(1-.02-.05)*pnorm((x-1)/.1),
    Quick=.02+(1-.02-.05)*(1-2^(-(x/1)^20)),
  ) %>% 
  pivot_longer(cols = !x) %>% 
  rename(Formulation=name)

pf_space_plot<-
  PFs %>% 
  ggplot()+
  geom_hline(yintercept = 1,color='white')+
  geom_line(aes(x=x,y=value,color=Formulation),linewidth=1,alpha=.8)+
  geom_vline(aes(xintercept = 0),linetype='dashed',color='grey')+

  theme_classic()+
  scale_color_manual(values = c('#009E73','#CC79A7'))+
  labs(
    title='Psychometric function formulation',
    x='Stimulus intensity (°C)',
    y='P(response="yes")',
    color=''
  )+
  scale_y_continuous(breaks = c(0,.5,1))

#Estimate and plot the two types of age-threshold relationships
ATR<-
  expand.grid(age=seq(20,80,1),idx=1:10^4) %>% 
  mutate(
    log_linear=exp(rnorm(n(),0+(age-20)*.05,.2))+rnorm(n()),
    linear=rnorm(n(),1+(age-20)*.1)
  ) %>% 
  pivot_longer(cols = !c(age,idx)) %>% 
  rename(Formulation=name) %>% 
  group_by(age,Formulation) %>% 
  summarise(m=median(value),lb=quantile(value,.025),ub=quantile(value,.975))
atr_space_plot<-
  ATR %>% 
  ggplot()+
  geom_ribbon(aes(x=age,ymin=lb,ymax=ub,fill = Formulation),alpha=.5)+
  geom_line(aes(x=age,y=m,color = Formulation))+
  theme_classic()+
  scale_color_manual(values = c('#4477AA','#EE6677'),labels=c('Linear','Log-linear'))+
  scale_fill_manual(values = c('#4477AA','#EE6677'),labels=c('Linear','Log-linear'))+
  labs(
    title='Effect of age formulation',
    x='Age (years)',
    y='Threshold (°C)',
    color='',
    fill=''
  )

#Combine and save
combined_model_space <- ggarrange(pf_space_plot,atr_space_plot, ncol = 1,heights = c(1,1))

ggsave(
  'plots/model_space.png',
  plot = combined_model_space,
  width = 10,
  height = 9,
  dpi = 300,
  units = 'cm'
)

#### Load fitted models ####
fit_a<-readRDS("results/fits/ageing/Gaussian_unconstrained_threshold_ageing_GLN.rds")
fit_an<-readRDS("results/fits/ageing_and_neuropathy/Gaussian_unconstrained_threshold_ageing_and_neuropathy_GLN_NI_H.rds")

#### Effect of age ####
#Extract age effect coefficients from the ageing only fit
age_effect<-
  fit_a$draws(variables = paste0('mu[',11:18,']'),format = 'df') %>%
  rename(
    CD_T=`mu[11]`,
    WD_T=`mu[12]`,
    CP_T=`mu[13]`,
    HP_T=`mu[14]`,
    CD_S=`mu[15]`,
    WD_S=`mu[16]`,
    CP_S=`mu[17]`,
    HP_S=`mu[18]`
  ) %>%
  pivot_longer(
    cols = c(CD_T, WD_T, CP_T, HP_T, CD_S, WD_S, CP_S, HP_S),
    names_to = c("task", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(task=factor(task,c('CD','WD','CP','HP'))) %>%
  select(task,parameter,value)

#Summarize posterior
P_A_T<-
  age_effect %>%
  filter(parameter=='T') %>%
  group_by(task) %>%
  summarise(
    es=sprintf("κ = %.2f [%.2f;%.2f]",mean(value),quantile(value,0.05),quantile(value,0.95)),
    p=sprintf("P(κ≤0) = %.3f",mean(value <= 0))
  )

P_A_S<-
  age_effect %>%
  filter(parameter=='S') %>%
  group_by(task) %>%
  summarise(
    es=sprintf("κ = %.2f [%.2f;%.2f]",mean(value),quantile(value,0.05),quantile(value,0.95)),
    p=sprintf("P(κ≥0) = %.3f",mean(value >= 0))
  )
#Sample from priors and compute Bayes factors
M<-8000
age_effect_priors<-
  tibble(
    CD_T=rnorm(M,0,0.54/60),
    WD_T=rnorm(M,0,1.90/60),
    CP_T=rnorm(M,0,7.39/60),
    HP_T=rnorm(M,0,2.67/60),
    CD_S=rnorm(M,0,1.11/60),
    WD_S=rnorm(M,0,0.73/60),
    CP_S=rnorm(M,0,0.50/60),
    HP_S=rnorm(M,0,0.79/60),
  ) %>%
  pivot_longer(
    cols = c(CD_T, WD_T, CP_T, HP_T, CD_S, WD_S, CP_S, HP_S),
    names_to = c("task", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(
    task=factor(task,c('CD','WD','CP','HP')),
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

#Estimate and summarize PFs
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
    g=`mu[9]`,
    l=`mu[10]`,
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
    guess=inv_logit(g),
    lapse=inv_logit(l),
    cdt=cd_ti+cd_ts*ra,
    wdt=wd_ti+wd_ts*ra,
    cpt=cp_ti+cp_ts*ra,
    hpt=hp_ti+hp_ts*ra,
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
    lb_60=quantile(p,0.20),
    ub_60=quantile(p,0.80),
    lb_80=quantile(p,0.10),
    ub_80=quantile(p,0.90),
    lb_90=quantile(p,0.05),
    ub_90=quantile(p,0.95),
    lb_95=quantile(p,0.025),
    ub_95=quantile(p,0.975)
  ) %>%
  pivot_longer(
    cols = matches("(lb|ub)_\\d+"),
    names_to = c("bound", "ci"),
    names_pattern = "(lb|ub)_(\\d+)",
    values_to = "value") %>%
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  )

pf_ageing$x[pf_ageing$task=='CD']<-30-pf_ageing$x[pf_ageing$task=='CD']
pf_ageing$x[pf_ageing$task=='WD']<-30+pf_ageing$x[pf_ageing$task=='WD']
pf_ageing$x[pf_ageing$task=='CP']<-30-pf_ageing$x[pf_ageing$task=='CP']
pf_ageing$x[pf_ageing$task=='HP']<-30+pf_ageing$x[pf_ageing$task=='HP']

#Make effect of age on threshold plot
ATP <- AE %>%
  filter(parameter == 'T') %>%
  ggplot(aes(x = value, color = task, linetype = distribution)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_density() +
  geom_point(
    data = dens_A %>% filter(parameter == "T"),
    aes(x = 0, y = density, shape = distribution,color=task),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = BF_A %>% filter(parameter == "T"),
    aes(x = 0.05, y = 48,label = BF01),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  geom_text(
    data = P_A_T,
    aes(x = 0.05, y = 54, label = p),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  geom_text(
    data = P_A_T,
    aes(x = 0.05, y = 60, label = es),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  scale_color_manual(
    values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('CD', 'WD', 'CP', 'HP')
  ) +
  guides(color='none',linetype='none',shape='none')+
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c(1,3)) +
  labs(
    title = 'Effect on the threshold',
    x = '',
    y = 'Density',
    color = 'Modality'
  ) +
  facet_wrap(~task, ncol = 4) +
  theme_classic() +
  theme(strip.text.x = element_blank()) +
  coord_cartesian(xlim = c(-.05, .3))

#Make effect of age on slope plot
ASP <- AE %>%
  filter(parameter == 'S') %>%
  ggplot(aes(x = value, color = task, linetype = distribution)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_density() +
  geom_point(
    data = dens_A %>% filter(parameter == "S"),
    aes(x = 0, y = density, shape = distribution,color=task),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = BF_A %>% filter(parameter == "S"),
    aes(x = -0.06,y=80,label = BF01),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  geom_text(
    data = P_A_S,
    aes(x = -0.06,y=90, label = p),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  geom_text(
    data = P_A_S,
    aes(x = -0.06, y = 100,label =es),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  scale_color_manual(
    values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('CD', 'WD', 'CP', 'HP')
  ) +
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c(1,3)) +
  labs(
    title = 'Effect on the slope (log-space)',
    x = '',
    y = 'Density',
    color = 'Modality: ',
    linetype='Distribution: ',
    shape='Distribution: '
  ) +
  facet_wrap(~task, ncol = 4) +
  theme_classic() +
  theme(strip.text.x = element_blank(),legend.position = 'bottom') +
  coord_cartesian(xlim = c(-.06,.03))

#Make group mean PF plot
PFAP<-pf_ageing %>%
  ggplot()+
  geom_vline(xintercept = 30, linetype = 'dashed', color = 'grey') +
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub,alpha=ci,fill=task))+
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
  scale_alpha_manual(
    values = c(.4,.3,.2,.1),
    labels=c('60%','80%','90%','95%'))+
  guides(color='none',fill='none')+
  xlab('Temperature (°C)')+
  ylab('P(response="yes")')+
  labs(
    fill='Modality',
    color='Modality',
    alpha='CI',
    title='Group mean psychometric functions at age 20 and 80'
  )+
  theme(
    strip.text.x = element_blank(),
    strip.background.x = element_blank(),
    legend.position = 'bottom'
  )+
  scale_y_continuous(breaks = c(0,.5,1))

#Combine and save
combinedA <- ggarrange(ATP, ASP, PFAP, ncol = 1,heights = c(1,1.2,2))

ggsave(
  'plots/ageing_effect_plot.png',
  plot = combinedA,
  width = 18,
  height = 24,
  dpi = 300,
  units = 'cm'
  )

#### Effect of neuropathy ####
#Extract neuropathy effect coefficients from the ageing and neuropathy fit
neuropathy_effect<-
  fit_an$draws(variables = paste0('mu[',11:18,']'),format = 'df') %>% 
  rename(
    CD_T=`mu[11]`,
    WD_T=`mu[12]`,
    CP_T=`mu[13]`,
    HP_T=`mu[14]`,
    CD_S=`mu[15]`,
    WD_S=`mu[16]`,
    CP_S=`mu[17]`,
    HP_S=`mu[18]`
  ) %>% 
  pivot_longer(
    cols = c(CD_T, WD_T, CP_T, HP_T, CD_S, WD_S, CP_S, HP_S),
    names_to = c("task", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>% 
  mutate(task=factor(task,c('CD','WD','CP','HP')))

#Summarize posterior
P_N_T<-
  neuropathy_effect %>%
  filter(parameter=='T') %>% 
  group_by(task) %>% 
  summarise(
    es=sprintf("κ = %.2f [%.2f;%.2f]",mean(value),quantile(value,0.05),quantile(value,0.95)),
    p=sprintf("P(κ≤0) = %.3f",mean(value <= 0))
  )
P_N_S<-
  neuropathy_effect %>%
  filter(parameter=='S') %>% 
  group_by(task) %>% 
  summarise(
    es=sprintf("κ = %.2f [%.2f;%.2f]",mean(value),quantile(value,0.05),quantile(value,0.95)),
    p=sprintf("P(κ≥0) = %.3f",mean(value >= 0))
  )
#Sample from priors and compute Bayes factors
M<-8000
neuropathy_effect_priors<-
  tibble(
    CD_T=rnorm(M,0,0.54),
    WD_T=rnorm(M,0,1.90),
    CP_T=rnorm(M,0,7.39),
    HP_T=rnorm(M,0,2.67),
    CD_S=rnorm(M,0,1.11),
    WD_S=rnorm(M,0,0.73),
    CP_S=rnorm(M,0,0.50),
    HP_S=rnorm(M,0,0.79),
    .iteration=1:M
  ) %>% 
  pivot_longer(
    cols = c(CD_T, WD_T, CP_T, HP_T, CD_S, WD_S, CP_S, HP_S),
    names_to = c("task", "parameter"),
    names_sep = "_",
    values_to = "value"
  ) %>% 
  mutate(
    task=factor(task,c('CD','WD','CP','HP')),
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
    g=`mu[9]`,
    l=`mu[10]`,
    cd_ts=`mu[11]`,
    wd_ts=`mu[12]`,
    cp_ts=`mu[13]`,
    hp_ts=`mu[14]`,
    cd_ss=`mu[15]`,
    wd_ss=`mu[16]`,
    cp_ss=`mu[17]`,
    hp_ss=`mu[18]`,
    cd_ta=`mu[19]`,
    wd_ta=`mu[20]`,
    cp_ta=`mu[21]`,
    hp_ta=`mu[22]`,
    cd_sa=`mu[23]`,
    wd_sa=`mu[24]`,
    cp_sa=`mu[25]`,
    hp_sa=`mu[26]`
  ) %>% 
  full_join(predictors) %>% 
  mutate(
    ra=age-20,
    guess=inv_logit(g),
    lapse=inv_logit(l),
    cdt=cd_ti+cd_ta*ra+cd_ts*status,
    wdt=wd_ti+wd_ta*ra+wd_ts*status,
    cpt=cp_ti+cp_ta*ra+cp_ts*status,
    hpt=hp_ti+hp_ta*ra+hp_ts*status,
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
    lb_60=quantile(p,0.20),    
    ub_60=quantile(p,0.80),
    lb_80=quantile(p,0.10),    
    ub_80=quantile(p,0.90),
    lb_90=quantile(p,0.05),    
    ub_90=quantile(p,0.95),
    lb_95=quantile(p,0.025),
    ub_95=quantile(p,0.975)
  ) %>% 
  pivot_longer(
    cols = matches("(lb|ub)_\\d+"),
    names_to = c("bound", "ci"), 
    names_pattern = "(lb|ub)_(\\d+)", 
    values_to = "value") %>% 
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  ) %>% 
  mutate(status_st=if_else(status==0, "Control", "Patient"))

pf_neuropathy$x[pf_neuropathy$task=='CD']<-30-pf_neuropathy$x[pf_neuropathy$task=='CD']
pf_neuropathy$x[pf_neuropathy$task=='WD']<-30+pf_neuropathy$x[pf_neuropathy$task=='WD']
pf_neuropathy$x[pf_neuropathy$task=='CP']<-30-pf_neuropathy$x[pf_neuropathy$task=='CP']
pf_neuropathy$x[pf_neuropathy$task=='HP']<-30+pf_neuropathy$x[pf_neuropathy$task=='HP']

#Make effect of neuropathy on threshold plot
NTP <- NE %>%
  filter(parameter == 'T') %>%
  ggplot(aes(x = value, color = task, linetype = distribution)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_density() +
  geom_point(
    data = dens_N %>% filter(parameter == "T"),
    aes(x = 0, y = density, shape = distribution,color=task),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = BF_N %>% filter(parameter == "T"),
    aes(x=1.5,y=.7,label = BF01),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  geom_text(
    data = P_N_T,
    aes(x=1.5,y=.8, label = p),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  geom_text(
    data = P_N_T,
    aes(x = 1.5, y = .9,label =  es),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  scale_color_manual(
    values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('CD', 'WD', 'CP', 'HP')
  ) +
  guides(color='none',linetype='none',shape='none')+
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c(1,3)) +
  labs(
    title = 'Effect on the threshold',
    x = '',
    y = 'Density',
    color = 'Modality'
  ) +
  facet_wrap(~task, ncol = 4) +
  theme_classic() +
  theme(strip.text.x = element_blank()) +
  coord_cartesian(xlim = c(-5, 11))

#Make effect of neuropathy on slope plot
NSP <- NE %>%
  filter(parameter == 'S') %>%
  ggplot(aes(x = value, color = task, linetype = distribution)) +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') +
  geom_density() +
  geom_point(
    data = dens_N %>% filter(parameter == "S"),
    aes(x = 0, y = density, shape = distribution,color=task),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = BF_N %>% filter(parameter == "S"),
    aes( x=-2,y=2.1,label = BF01),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  geom_text(
    data = P_N_S,
    aes( x=-2,y=2.3, label = p),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  geom_text(
    data = P_N_S,
    aes(x = -2, y = 2.5,label = es),
    inherit.aes = FALSE,
    hjust = 'left', size = 2.5,
    parse = F
  ) +
  scale_color_manual(
    values = c('#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('CD', 'WD', 'CP', 'HP')
  ) +
  scale_shape_manual(values = c(16, 1)) +
  scale_linetype_manual(values = c(1,3)) +
  labs(
    title = 'Effect on the slope (log-space)',
    x = '',
    y = 'Density',
    color = 'Modality: ',
    linetype='Distribution: ',
    shape='Distribution: '
  ) +
  facet_wrap(~task, ncol = 4) +
  theme_classic() +
  theme(strip.text.x = element_blank(),legend.position = 'bottom') +
  coord_cartesian(xlim = c(-2,1))

#Make group mean PF plot
PFNP<-pf_neuropathy %>%
  ggplot()+
  geom_vline(xintercept = 30, linetype = 'dashed', color = 'grey') +
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub,alpha=ci,fill=task))+
  geom_line(aes(x=x,y=m,color=task))+
  theme_classic()+
  ggh4x::facet_grid2(
    rows = vars(status_st),
    cols = vars(task),
    scales = "free",
    axes = "x"
  ) +  
  scale_color_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  scale_fill_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  scale_alpha_manual(
    values = c(.4,.3,.2,.1),
    labels=c('60%','80%','90%','95%'))+
  guides(color='none',fill='none')+
  xlab('Temperature (°C)')+
  ylab('P(response="yes")')+
  labs(
    fill='Modality',
    color='Modality',
    alpha='CI',
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
print(PFNP)

#Combine and save
combinedN <- ggarrange(NTP, NSP, PFNP, ncol = 1,heights = c(1,1.2,2))

ggsave(
  'plots/neuropathy_effect_plot.png',
  plot = combinedN,
  width = 18,
  height = 24,
  dpi = 300,
  units = 'cm'
  )

#### Discrimination ####
#load discrimination analysis results
s_and_s<-readRDS('results/S_and_S_FUP.rds')
AUC<-readRDS('results/AUC_FUP.rds')

#Reformat and summarize
temp<-
  s_and_s %>% 
  mutate(
    model=factor(model,labels = c('CDT','CDS','CDT&CDS','WDT','WDS','WDT&WDS','CPT','CPS','CPT&CPS','HPT','HPS','HPT&HPS','All')),
    specificity=round(specificity,2)
  ) %>% 
  group_by(model,specificity) %>%
  summarise(
    m=median(sensitivity),
    lb_60=quantile(sensitivity,0.20),
    ub_60=quantile(sensitivity,0.80),
    lb_80=quantile(sensitivity,0.10),
    ub_80=quantile(sensitivity,0.90),
    lb_90=quantile(sensitivity,0.05),
    ub_90=quantile(sensitivity,0.95)
  ) %>%
  pivot_longer(
    cols = matches("(lb|ub)_\\d+"),
    names_to = c("bound", "ci"),
    names_pattern = "(lb|ub)_(\\d+)",
    values_to = "value") %>%
  pivot_wider(
    names_from = "bound",
    values_from = "value"
  )

temp$task[temp$model%in%c('CDT','CDS','CDT&CDS')]<-'CD'
temp$task[temp$model%in%c('WDT','WDS','WDT&WDS')]<-'WD'
temp$task[temp$model%in%c('CPT','CPS','CPT&CPS')]<-'CP'
temp$task[temp$model%in%c('HPT','HPS','HPT&HPS')]<-'HP'
temp$task[temp$model%in%c('All')]<-'All'
temp$param[temp$model%in%c('CDT','WDT','CPT','HPT')]<-'Threshold'
temp$param[temp$model%in%c('CDS','WDS','CPS','HPS')]<-'Slope'
temp$param[temp$model%in%c('CDT&CDS','WDT&WDS','CPT&CPS','HPT&HPS')]<-'Both'
temp$param[temp$model%in%c('All')]<-'All'

temp$param<-factor(temp$param,levels=c('Threshold','Slope','Both','All'))
temp$task<-factor(temp$task,levels=c('CD','WD','CP','HP'))

auroc<-
  AUC %>% 
  group_by(task,param,model) %>% 
  summarise(m=round(median(AUC),2),
            lb=round(quantile(AUC,0.025),2),
            ub=round(quantile(AUC,.975),2),
            p=mean(AUC<=0.5)
  ) %>% 
  mutate(
    AUC=sprintf('AUC = %.2f [%.2f;%.2f]',m,lb,ub),
    p=sprintf('P(AUC≤0.5) = %.3f',p)
    )
auroc_differences<-
  AUC %>% 
  filter(!is.na(task)) %>% 
  select(!model) %>% 
  pivot_wider(names_from = param,values_from = AUC) %>% 
  mutate(
    t_minus_s=Threshold-Slope,
    t_minus_b=Threshold-Both,
    b_minus_s=Both-Slope
  ) %>% 
  pivot_longer(cols = c(t_minus_s,t_minus_b,b_minus_s)) %>% 
  group_by(task,name) %>% 
  summarise(m=round(median(value),2),
            lb=round(quantile(value,0.025),2),
            ub=round(quantile(value,.975),2),
            p=mean(value<=0)
  )
#Make plot and save it
temp %>% 
  filter(param!='All') %>% 
  ggplot()+
  geom_abline()+
  geom_ribbon(aes(x=1-specificity,ymin=lb,ymax=ub,alpha=ci,fill=task))+
  geom_line(aes(x=1-specificity,y=m,color=task))+
  geom_text(
    data = filter(auroc,param!='All'),
    aes(x=1,y=.2,label = AUC),
    hjust='right',
    size=2
  )+
  geom_text(
    data = filter(auroc,param!='All'),
    aes(x=1,y=.1,label = p),
    hjust='right',
    size=2
  )+
  theme_classic()+
  facet_grid(cols = vars(task),rows=vars(param)) +   
  scale_color_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  scale_fill_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  scale_alpha_manual(
    values = c(.4,.3,.2),
    labels=c('60%','80%','90%'))+
  xlab('1-Specificty')+
  ylab('Sensitivity')+
  labs(
    fill='Modality',
    color='Modality',
    alpha='CI',
    title='Receiver Operating Characteristic Curves'
  )+
  theme(
    strip.text.x = element_blank(),
    strip.background.x = element_blank()
  )+
  scale_y_continuous(breaks = c(0,1))+
  scale_x_continuous(breaks = c(0,1))

ggsave(
  'plots/ROC_FUP.png',
  width = 18,
  height = 10,
  dpi = 300,
  units = 'cm'
)

#### Correlation ####
#Extract correlation matrix from ageing and neuropathy fit
draws_df <- fit_an$draws(variables = "omega", format = "draws_df")
correlation_draws_long <- draws_df %>%
  select(starts_with("omega[")) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(
    i = as.integer(str_extract(param, "(?<=\\[)\\d+")),
    j = as.integer(str_extract(param, "(?<=,)\\d+"))
  )

correlation_draws_filtered <- correlation_draws_long %>%
  filter(i > j,i<9,j<9)

correlation_draws_filtered$i<-factor(correlation_draws_filtered$i,labels=c('WDT','CPT','HPT','CDS','WDS','CPS','HPS'))
correlation_draws_filtered$j<-factor(correlation_draws_filtered$j,labels=c('CDT','WDT','CPT','HPT','CDS','WDS','CPS'))

#Summarize
rhos<-
  correlation_draws_filtered %>% 
  group_by(i,j) %>% 
  summarise(
    m=sprintf("%.2f",median(value)),
    lb=sprintf("%.2f",quantile(value,0.025)),
    ub=sprintf("%.2f",quantile(value,0.975))
  )

#Make plot and save
corr_plot<-
  correlation_draws_filtered %>% 
  ggplot() +
  geom_vline(aes(xintercept = -1),linetype='dashed',color='grey')+
  geom_vline(aes(xintercept = 0),linetype='dashed',color='grey')+
  geom_vline(aes(xintercept = 1),linetype='dashed',color='grey')+
  geom_density(aes(x = value),color = "black", fill = "white") +
  geom_text(
    data=rhos,
    aes(x=-0.5,y=5.2,label=paste0(m,'\n[',lb,';',ub,']')),
    size=2,
    fontface = "bold"
  )+
  facet_grid(i ~ j) +
  theme_classic() +
  labs(
    title = 'Partial correlations between parameters',
    x = "Correlation", 
    y = ""
  )+
  scale_x_continuous(breaks = -1:1)+
  scale_y_continuous(breaks = c())+
  coord_cartesian(ylim=c(0,6))

ggsave(
  'plots/corr_plot_mod.png',
  plot = corr_plot,
  width = 18,
  height = 18,
  dpi = 300,
  units = 'cm'
)

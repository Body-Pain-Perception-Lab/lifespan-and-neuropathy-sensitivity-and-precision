# Script to generate single subject postrior predictive plots
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)

rm(list=ls())

#### Load fitted models ####
fit_an<-readRDS("results/fits/Gaussian_constrained_threshold_ageing_and_neuropathy_NRS_NI.rds")

#Extract participant data
d <- read_csv("data/aggregated_results.csv") %>% 
  filter(
    recording_deviates_from_target==0,
    recording_deviates_from_mean==0,
    !is.na(response),
    !subject%in%c(3019), #exclusion based on saturation check plots
  ) %>% 
  mutate(recorded_intensity=round(recorded_intensity,1)) %>% 
  group_by(recorded_intensity,subject,condition) %>% 
  summarise(prop=mean(response),n=sum(response>-1)) %>% 
  arrange(subject)

d$task[d$condition==1]<-'CD'
d$task[d$condition==2]<-'WD'
d$task[d$condition==3]<-'CP'
d$task[d$condition==4]<-'HP'

participant<-unique(d$subject)

for(pdx in 1:length(unique(d$subject))){
  d$subject[d$subject==participant[pdx]]<-pdx
}
d<-d %>% rename(participant=subject,x=recorded_intensity)

#Extract participant parameters
predictors<-expand_grid(x=seq(0.1,30,.1),.draw=1:8000)

draws_an<-fit_an$draws(variables = c('alpha','beta','gamma','lambda'),format='df')%>%
  pivot_longer(!starts_with('.'), names_to = "param", values_to = "value")

draws_an_task <- draws_an %>%
  filter(str_detect(param, "alpha|beta")) %>%
  extract(param, into = c("param", "task", "participant"),
          regex = "([a-z]+)\\[([^,]+),([^\\]]+)\\]") %>%
  pivot_wider(names_from = param, values_from = value)

draws_an_nontask <- draws_an %>%
  filter(str_detect(param, "gamma|lambda")) %>%
  extract(param, into = c("param", "participant"),
          regex = "([a-z]+)\\[([^\\]]+)\\]") %>%
  pivot_wider(names_from = param, values_from = value)

# Join
draws_an <- draws_an_task %>%
  full_join(draws_an_nontask) %>% 
  mutate(participant= as.integer(participant)) %>%
  arrange(participant) %>% 
  full_join(predictors) %>% 
  mutate(p=gamma+(1-gamma-lambda)*pnorm(x,alpha,1/beta)) %>% 
  group_by(x,participant,task) %>% 
  summarise(
    m=median(p),
    lb=quantile(p,0.025),
    ub=quantile(p,0.975)
  )
draws_an$task[draws_an$task=='1']<-'CD'
draws_an$task[draws_an$task=='2']<-'WD'
draws_an$task[draws_an$task=='3']<-'CP'
draws_an$task[draws_an$task=='4']<-'HP'

# Plot CDT
CD_fig<-
  draws_an %>% 
  filter(task=='CD') %>% 
  ggplot()+
  geom_vline(aes(xintercept=0),'linetype'='dotted')+
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub),fill='#56B4E9',alpha=.5)+
  geom_line(aes(x=x,y=m),color='#56B4E9')+
  geom_point(data=d %>% filter(task=='CD'),aes(x=x,y=prop,color=n))+
  theme_minimal()+
  facet_wrap(participant~.,ncol = 9)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab('Temperature (°C)')+
  ylab('P(response="yes")')+
  labs(
    alpha='CI',
    title='Participant fits - cold detection'
  )+
  scale_color_continuous(trans = "log10",type = 'viridis')
ggsave(
  'plots/pp_cd.jpeg',
  plot = CD_fig,
  width = 18,
  height = 24,
  dpi = 300,
  units = 'cm'
)

CP_fig<-
  draws_an %>% 
  filter(task=='CP') %>% 
  ggplot()+
  geom_vline(aes(xintercept=0),'linetype'='dotted')+
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub),fill='#0072B2',alpha=.5)+
  geom_line(aes(x=x,y=m),color='#0072B2')+
  geom_point(data=d %>% filter(task=='CP'),aes(x=x,y=prop,color=n))+
  theme_minimal()+
  facet_wrap(participant~.,ncol = 9)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab('Temperature (°C)')+
  ylab('P(response="yes")')+
  labs(
    alpha='CI',
    title='Participant fits - cold pain'
  )+
  scale_color_continuous(trans = "log10",type = 'viridis')
ggsave(
  'plots/pp_cp.jpeg',
  plot = CP_fig,
  width = 18,
  height = 24,
  dpi = 300,
  units = 'cm'
)

WD_fig<-
  draws_an %>% 
  filter(task=='WD') %>% 
  ggplot()+
  geom_vline(aes(xintercept=0),'linetype'='dotted')+
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub),fill='#E69F00',alpha=.5)+
  geom_line(aes(x=x,y=m),color='#E69F00')+
  geom_point(data=d %>% filter(task=='WD'),aes(x=x,y=prop,color=n))+
  theme_minimal()+
  facet_wrap(participant~.,ncol = 9)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab('Temperature (°C)')+
  ylab('P(response="yes")')+
  labs(
    alpha='CI',
    title='Participant fits - warm detection'
  )+
  scale_color_continuous(trans = "log10",type = 'viridis')
ggsave(
  'plots/pp_wd.jpeg',
  plot = WD_fig,
  width = 18,
  height = 24,
  dpi = 300,
  units = 'cm'
)

HP_fig<-
  draws_an %>% 
  filter(task=='HP') %>% 
  ggplot()+
  geom_vline(aes(xintercept=0),'linetype'='dotted')+
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub),fill='#D55E00',alpha=.5)+
  geom_line(aes(x=x,y=m),color='#D55E00')+
  geom_point(data=d %>% filter(task=='HP'),aes(x=x,y=prop,color=n))+
  theme_minimal()+
  facet_wrap(participant~.,ncol = 9)+
  scale_y_continuous(breaks=c(0,0.5,1))+
  xlab('Temperature (°C)')+
  ylab('P(response="yes")')+
  labs(
    alpha='CI',
    title='Participant fits - heat pain'
  )+
  scale_color_continuous(trans = "log10",type = 'viridis')
ggsave(
  'plots/pp_hp.jpeg',
  plot = HP_fig,
  width = 18,
  height = 24,
  dpi = 300,
  units = 'cm'
)

# Script to generate figure 3
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(ggpubr)
rm(list=ls())

#### Model space illustration ####
#Estimate and plot the two types of PF
PFs<-
  tibble(x=seq(0.01,2,.01)) %>% 
  mutate(
    Gaussian=.02+(1-.02-.05)*pnorm((x-1)/.1),
  ) %>% 
  pivot_longer(cols = !x) %>% 
  rename(Formulation=name)

pf_plot<-
  PFs %>% 
  ggplot()+
  geom_hline(yintercept = 1,color='white')+
  geom_line(aes(x=x,y=value),linewidth=1,alpha=.8,color='#009E73')+
  geom_vline(aes(xintercept = 0),linetype='dashed',color='grey')+
  theme_classic()+
  labs(
    title='Psychometric function formulation',
    x='Stimulus intensity',
    y='P(response="yes")',
    color=''
  )+
  scale_y_continuous(breaks = c(0,.5,1))

ggsave(
  'plots/model.jpeg',
  plot = pf_plot,
  width = 10,
  height = 5,
  dpi = 300,
  units = 'cm'
)

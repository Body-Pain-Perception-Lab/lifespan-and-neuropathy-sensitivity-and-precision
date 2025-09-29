# Script to generate prior predictive plots from the weakly informative priors 
# and the priors informed by a previous dataset
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set up environement ####
library(tidyverse)

inv_logit<-function(x){
  y=1/(1+exp(-x))
  return(y)
}

#### Gaussian Unconstrained ####
M<-10^4

#Sample group parameters from the weakly informative prior
Initial_mu_1<-rnorm(M,5,5)
Initial_mu_2<-rnorm(M,5,5)
Initial_mu_3<-rnorm(M,15,10)
Initial_mu_4<-rnorm(M,15,10)
Initial_mu_5<-rnorm(M)
Initial_mu_6<-rnorm(M)
Initial_mu_7<-rnorm(M)
Initial_mu_8<-rnorm(M)
Initial_mu_9<-rnorm(M,-4,0.5)
Initial_mu_10<-rnorm(M,-4,0.5)

Initial_tau_1<-abs(rnorm(M,0,5))
Initial_tau_2<-abs(rnorm(M,0,5))
Initial_tau_3<-abs(rnorm(M,0,10))
Initial_tau_4<-abs(rnorm(M,0,10))
Initial_tau_5<-abs(rnorm(M))
Initial_tau_6<-abs(rnorm(M))
Initial_tau_7<-abs(rnorm(M))
Initial_tau_8<-abs(rnorm(M))
Initial_tau_9<-abs(rnorm(M,0,0.5))
Initial_tau_10<-abs(rnorm(M,0,0.5))

Initial_kappa1<-Initial_kappa2<-Initial_kappa3<-Initial_kappa4<-Initial_kappa5<-
  Initial_kappa6<-Initial_kappa7<-Initial_kappa8<-Initial_kappa9<-Initial_kappa10<-
  vector(,M)

for(m in 1:M){
  Initial_kappa1[m]<-rnorm(1,Initial_mu_1[m],Initial_tau_1[m])
  Initial_kappa2[m]<-rnorm(1,Initial_mu_2[m],Initial_tau_2[m])
  Initial_kappa3[m]<-rnorm(1,Initial_mu_3[m],Initial_tau_3[m])
  Initial_kappa4[m]<-rnorm(1,Initial_mu_4[m],Initial_tau_4[m])
  Initial_kappa5[m]<-rnorm(1,Initial_mu_5[m],Initial_tau_5[m])
  Initial_kappa6[m]<-rnorm(1,Initial_mu_6[m],Initial_tau_6[m])
  Initial_kappa7[m]<-rnorm(1,Initial_mu_7[m],Initial_tau_7[m])
  Initial_kappa8[m]<-rnorm(1,Initial_mu_8[m],Initial_tau_8[m])
  Initial_kappa9[m]<-rnorm(1,Initial_mu_9[m],Initial_tau_9[m])
  Initial_kappa10[m]<-rnorm(1,Initial_mu_10[m],Initial_tau_10[m])
}

#Sample group parameters from the informed prior
Final_mu_1<-rnorm(M,0.67,0.09*4)
Final_mu_2<-rnorm(M,4.33,0.31*4)
Final_mu_3<-rnorm(M,22.48,1.18*4)
Final_mu_4<-rnorm(M,15.11,0.42*4)
Final_mu_5<-rnorm(M,1.19,0.22*4)
Final_mu_6<-rnorm(M,-0.23,0.15*4)
Final_mu_7<-rnorm(M,-1.60,0.11*4)
Final_mu_8<-rnorm(M,-0.22,0.15*4)
Final_mu_9<-rnorm(M,-3.33,0.40)
Final_mu_10<-rnorm(M,-3.84,0.36)

Final_tau_1<-abs(rnorm(M,0.54,0.09*4))
Final_tau_2<-abs(rnorm(M,1.90,0.24*4))
Final_tau_3<-abs(rnorm(M,7.39,0.83*4))
Final_tau_4<-abs(rnorm(M,2.67,0.30*4))
Final_tau_5<-abs(rnorm(M,1.11,0.23*4))
Final_tau_6<-abs(rnorm(M,0.73,0.13*4))
Final_tau_7<-abs(rnorm(M,0.50,0.13*4))
Final_tau_8<-abs(rnorm(M,0.79,0.14*4))
Final_tau_9<-abs(rnorm(M,0.66,0.35))
Final_tau_10<-abs(rnorm(M,1.13,0.45))

Final_kappa1<-Final_kappa2<-Final_kappa3<-Final_kappa4<-Final_kappa5<-
  Final_kappa6<-Final_kappa7<-Final_kappa8<-Final_kappa9<-Final_kappa10<-
  vector(,M)

for(m in 1:M){
  Final_kappa1[m]<-rnorm(1,Final_mu_1[m],Final_tau_1[m])
  Final_kappa2[m]<-rnorm(1,Final_mu_2[m],Final_tau_2[m])
  Final_kappa3[m]<-rnorm(1,Final_mu_3[m],Final_tau_3[m])
  Final_kappa4[m]<-rnorm(1,Final_mu_4[m],Final_tau_4[m])
  Final_kappa5[m]<-rnorm(1,Final_mu_5[m],Final_tau_5[m])
  Final_kappa6[m]<-rnorm(1,Final_mu_6[m],Final_tau_6[m])
  Final_kappa7[m]<-rnorm(1,Final_mu_7[m],Final_tau_7[m])
  Final_kappa8[m]<-rnorm(1,Final_mu_8[m],Final_tau_8[m])
  Final_kappa9[m]<-rnorm(1,Final_mu_9[m],Final_tau_9[m])
  Final_kappa10[m]<-rnorm(1,Final_mu_10[m],Final_tau_10[m])
}

#Combine and reformat draws
g<-expand.grid(x=seq(0,35,.2),idx=1:M)
kappas<-tibble(
  idx=1:M,
  Initial_kappa1=Initial_kappa1,
  Initial_kappa2=Initial_kappa2,
  Initial_kappa3=Initial_kappa3,
  Initial_kappa4=Initial_kappa4,
  Initial_kappa5=Initial_kappa5,
  Initial_kappa6=Initial_kappa6,
  Initial_kappa7=Initial_kappa7,
  Initial_kappa8=Initial_kappa8,
  Initial_kappa9=Initial_kappa9,
  Initial_kappa10=Initial_kappa10,
  Final_kappa1=Final_kappa1,
  Final_kappa2=Final_kappa2,
  Final_kappa3=Final_kappa3,
  Final_kappa4=Final_kappa4,
  Final_kappa5=Final_kappa5,
  Final_kappa6=Final_kappa6,
  Final_kappa7=Final_kappa7,
  Final_kappa8=Final_kappa8,
  Final_kappa9=Final_kappa9,
  Final_kappa10=Final_kappa10
)

#Estimate and summarize PFs
PFS<-
  full_join(g,kappas)%>%
  pivot_longer(
    cols = !c(x,idx),
    names_to = c('stage','kappa'),
    names_pattern = "(.*)_(.*)"
  ) %>% 
  pivot_wider(
    names_from = kappa
  ) %>% 
  mutate(
    CD=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*pnorm(x,(kappa1),1/exp(kappa5)),    
    WD=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*pnorm(x,(kappa2),1/exp(kappa6)),    
    CP=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*pnorm(x,(kappa3),1/exp(kappa7)),    
    HP=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*pnorm(x,(kappa4),1/exp(kappa8))
  ) %>% 
  pivot_longer(
    cols = c(CD,CP,WD,HP)
  ) %>% 
  group_by(name,x,stage) %>% 
  summarise(
    m=median(value),
    lb_60=quantile(value,0.20),
    ub_60=quantile(value,0.80),    
    lb_80=quantile(value,0.10),
    ub_80=quantile(value,0.90),    
    lb_90=quantile(value,0.05),
    ub_90=quantile(value,0.95),
    lb_95=quantile(value,0.025),
    ub_95=quantile(value,0.975)
  ) %>% 
  pivot_longer(
    cols =!c(m,x,name,stage),
    names_to = c('bound','CI'),
    names_pattern = "(.*)_(.*)"
  ) %>% 
  pivot_wider(
    names_from = bound,
    values_from = value
  ) %>% 
  filter((name=='CD'&x<15)|(name=='WD'&x<15)|name=='CP'|(name=='HP'&x<20))

#Make and save the plot
PFS$name<-factor(PFS$name,c('CD','WD','CP','HP'))
PFS$stage<-factor(PFS$stage,c('Initial','Final'))
PFS$x[PFS$name=='CD']<-30-PFS$x[PFS$name=='CD']
PFS$x[PFS$name=='WD']<-30+PFS$x[PFS$name=='WD']
PFS$x[PFS$name=='CP']<-30-PFS$x[PFS$name=='CP']
PFS$x[PFS$name=='HP']<-30+PFS$x[PFS$name=='HP']

PFS %>% 
  ggplot()+
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub,alpha=CI,fill=name))+
  geom_line(aes(x=x,y=m,color=name))+
  facet_grid(stage~name,scales = 'free')+
  labs(
    title='Gaussian - unconstrained threshold',
    subtitle='Prior predictive distribution',
    x='Temperature (°C)',
    y='P(response="yes")',
    fill='Modality',
    color='Modality',
    alpha='Quantiles',
  )+
  scale_alpha_manual(labels=c('60%','80%','90%','95%'),values=c(.4,.3,.2,.1))+
  scale_color_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  scale_fill_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  theme_classic() +
  theme(strip.text.x = element_blank())

ggsave('GUpriors.png',units = 'cm',width = 16,height = 10)

#### Gaussian Constrained ####

M<-10^4

#Sample group parameters from the weakly informative prior
Initial_mu_1<-rnorm(M)
Initial_mu_2<-rnorm(M)
Initial_mu_3<-rnorm(M,2.5,0.5)
Initial_mu_4<-rnorm(M,2.5,0.5)
Initial_mu_5<-rnorm(M)
Initial_mu_6<-rnorm(M)
Initial_mu_7<-rnorm(M)
Initial_mu_8<-rnorm(M)
Initial_mu_9<-rnorm(M,-4,0.5)
Initial_mu_10<-rnorm(M,-4,0.5)

Initial_tau_1<-abs(rnorm(M))
Initial_tau_2<-abs(rnorm(M))
Initial_tau_3<-abs(rnorm(M,0,0.5))
Initial_tau_4<-abs(rnorm(M,0,0.5))
Initial_tau_5<-abs(rnorm(M))
Initial_tau_6<-abs(rnorm(M))
Initial_tau_7<-abs(rnorm(M))
Initial_tau_8<-abs(rnorm(M))
Initial_tau_9<-abs(rnorm(M,0,0.5))
Initial_tau_10<-abs(rnorm(M,0,0.5))

Initial_kappa1<-Initial_kappa2<-Initial_kappa3<-Initial_kappa4<-Initial_kappa5<-
  Initial_kappa6<-Initial_kappa7<-Initial_kappa8<-Initial_kappa9<-Initial_kappa10<-
  vector(,M)

for(m in 1:M){
  Initial_kappa1[m]<-rnorm(1,Initial_mu_1[m],Initial_tau_1[m])
  Initial_kappa2[m]<-rnorm(1,Initial_mu_2[m],Initial_tau_2[m])
  Initial_kappa3[m]<-rnorm(1,Initial_mu_3[m],Initial_tau_3[m])
  Initial_kappa4[m]<-rnorm(1,Initial_mu_4[m],Initial_tau_4[m])
  Initial_kappa5[m]<-rnorm(1,Initial_mu_5[m],Initial_tau_5[m])
  Initial_kappa6[m]<-rnorm(1,Initial_mu_6[m],Initial_tau_6[m])
  Initial_kappa7[m]<-rnorm(1,Initial_mu_7[m],Initial_tau_7[m])
  Initial_kappa8[m]<-rnorm(1,Initial_mu_8[m],Initial_tau_8[m])
  Initial_kappa9[m]<-rnorm(1,Initial_mu_9[m],Initial_tau_9[m])
  Initial_kappa10[m]<-rnorm(1,Initial_mu_10[m],Initial_tau_10[m])
}

#Sample group parameters from the informed prior
Final_mu_1<-rnorm(M,-0.52,0.15*4)
Final_mu_2<-rnorm(M,1.35,0.06*4)
Final_mu_3<-rnorm(M,3.06,0.08*4)
Final_mu_4<-rnorm(M,2.71,0.03*4)
Final_mu_5<-rnorm(M,1.20,0.20*4)
Final_mu_6<-rnorm(M,-0.24,0.15*4)
Final_mu_7<-rnorm(M,-1.64,0.12*4)
Final_mu_8<-rnorm(M,-0.30,0.16*4)
Final_mu_9<-rnorm(M,-3.32,0.41)
Final_mu_10<-rnorm(M,-3.74,0.35)

Final_tau_1<-abs(rnorm(M,0.92,0.13*4))
Final_tau_2<-abs(rnorm(M,0.53,0.07*4))
Final_tau_3<-abs(rnorm(M,0.38,0.04*4))
Final_tau_4<-abs(rnorm(M,0.17,0.04*4))
Final_tau_5<-abs(rnorm(M,1.09,0.19*4))
Final_tau_6<-abs(rnorm(M,0.75,0.14*4))
Final_tau_7<-abs(rnorm(M,0.53,0.11*4))
Final_tau_8<-abs(rnorm(M,0.85,0.14*4))
Final_tau_9<-abs(rnorm(M,0.63,0.36))
Final_tau_10<-abs(rnorm(M,0.55,0.38))

Final_kappa1<-Final_kappa2<-Final_kappa3<-Final_kappa4<-Final_kappa5<-
  Final_kappa6<-Final_kappa7<-Final_kappa8<-Final_kappa9<-Final_kappa10<-
  vector(,M)

for(m in 1:M){
  Final_kappa1[m]<-rnorm(1,Final_mu_1[m],Final_tau_1[m])
  Final_kappa2[m]<-rnorm(1,Final_mu_2[m],Final_tau_2[m])
  Final_kappa3[m]<-rnorm(1,Final_mu_3[m],Final_tau_3[m])
  Final_kappa4[m]<-rnorm(1,Final_mu_4[m],Final_tau_4[m])
  Final_kappa5[m]<-rnorm(1,Final_mu_5[m],Final_tau_5[m])
  Final_kappa6[m]<-rnorm(1,Final_mu_6[m],Final_tau_6[m])
  Final_kappa7[m]<-rnorm(1,Final_mu_7[m],Final_tau_7[m])
  Final_kappa8[m]<-rnorm(1,Final_mu_8[m],Final_tau_8[m])
  Final_kappa9[m]<-rnorm(1,Final_mu_9[m],Final_tau_9[m])
  Final_kappa10[m]<-rnorm(1,Final_mu_10[m],Final_tau_10[m])
}

#Combine and reformat draws
g<-expand.grid(x=seq(0,35,.2),idx=1:M)
kappas<-tibble(
  idx=1:M,
  Initial_kappa1=Initial_kappa1,
  Initial_kappa2=Initial_kappa2,
  Initial_kappa3=Initial_kappa3,
  Initial_kappa4=Initial_kappa4,
  Initial_kappa5=Initial_kappa5,
  Initial_kappa6=Initial_kappa6,
  Initial_kappa7=Initial_kappa7,
  Initial_kappa8=Initial_kappa8,
  Initial_kappa9=Initial_kappa9,
  Initial_kappa10=Initial_kappa10,
  Final_kappa1=Final_kappa1,
  Final_kappa2=Final_kappa2,
  Final_kappa3=Final_kappa3,
  Final_kappa4=Final_kappa4,
  Final_kappa5=Final_kappa5,
  Final_kappa6=Final_kappa6,
  Final_kappa7=Final_kappa7,
  Final_kappa8=Final_kappa8,
  Final_kappa9=Final_kappa9,
  Final_kappa10=Final_kappa10
)

#Estimate and summarize PFs
PFS<-
  full_join(g,kappas)%>%
  pivot_longer(
    cols = !c(x,idx),
    names_to = c('stage','kappa'),
    names_pattern = "(.*)_(.*)"
  ) %>% 
  pivot_wider(
    names_from = kappa
  ) %>% 
  mutate(
    CD=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*pnorm(x,exp(kappa1),1/exp(kappa5)),    
    WD=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*pnorm(x,exp(kappa2),1/exp(kappa6)),    
    CP=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*pnorm(x,exp(kappa3),1/exp(kappa7)),    
    HP=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*pnorm(x,exp(kappa4),1/exp(kappa8))
  ) %>% 
  pivot_longer(
    cols = c(CD,CP,WD,HP)
  ) %>% 
  group_by(name,x,stage) %>% 
  summarise(
    m=median(value),
    lb_60=quantile(value,0.20),
    ub_60=quantile(value,0.80),    
    lb_80=quantile(value,0.10),
    ub_80=quantile(value,0.90),    
    lb_90=quantile(value,0.05),
    ub_90=quantile(value,0.95),
    lb_95=quantile(value,0.025),
    ub_95=quantile(value,0.975)
  ) %>% 
  pivot_longer(
    cols =!c(m,x,name,stage),
    names_to = c('bound','CI'),
    names_pattern = "(.*)_(.*)"
  ) %>% 
  pivot_wider(
    names_from = bound,
    values_from = value
  ) %>% 
  filter((name=='CD'&x<15)|(name=='WD'&x<15)|name=='CP'|(name=='HP'&x<20))

#Make and save the plot
PFS$name<-factor(PFS$name,c('CD','WD','CP','HP'))
PFS$stage<-factor(PFS$stage,c('Initial','Final'))
PFS$x[PFS$name=='CD']<-30-PFS$x[PFS$name=='CD']
PFS$x[PFS$name=='WD']<-30+PFS$x[PFS$name=='WD']
PFS$x[PFS$name=='CP']<-30-PFS$x[PFS$name=='CP']
PFS$x[PFS$name=='HP']<-30+PFS$x[PFS$name=='HP']

PFS %>% 
  ggplot()+
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub,alpha=CI,fill=name))+
  geom_line(aes(x=x,y=m,color=name))+
  facet_grid(stage~name,scales = 'free')+
  labs(
    title='Gaussian - constrained threshold',
    subtitle='Prior predictive distribution',
    x='Temperature (°C)',
    y='P(response="yes")',
    fill='Modality',
    color='Modality',
    alpha='Quantiles',
  )+
  scale_alpha_manual(labels=c('60%','80%','90%','95%'),values=c(.4,.3,.2,.1))+
  scale_color_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  scale_fill_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  theme_classic() +
  theme(strip.text.x = element_blank())

ggsave('GCpriors.png',units = 'cm',width = 16,height = 10)

#### Quick Constrained ####

M<-10^4

#Sample group parameters from the weakly informative prior
Initial_mu_1<-rnorm(M)
Initial_mu_2<-rnorm(M)
Initial_mu_3<-rnorm(M,2.5,0.5)
Initial_mu_4<-rnorm(M,2.5,0.5)
Initial_mu_5<-rnorm(M,1.5,1.5)
Initial_mu_6<-rnorm(M,1.5,1.5)
Initial_mu_7<-rnorm(M,2.5,1.25)
Initial_mu_8<-rnorm(M,2.5,1.25)
Initial_mu_9<-rnorm(M,-4,0.5)
Initial_mu_10<-rnorm(M,-4,0.5)

Initial_tau_1<-abs(rnorm(M))
Initial_tau_2<-abs(rnorm(M))
Initial_tau_3<-abs(rnorm(M,0,0.5))
Initial_tau_4<-abs(rnorm(M,0,0.5))
Initial_tau_5<-abs(rnorm(M,0,1.5))
Initial_tau_6<-abs(rnorm(M,0,1.5))
Initial_tau_7<-abs(rnorm(M,0,1.25))
Initial_tau_8<-abs(rnorm(M,0,1.25))
Initial_tau_9<-abs(rnorm(M,0,0.5))
Initial_tau_10<-abs(rnorm(M,0,0.5))

Initial_kappa1<-Initial_kappa2<-Initial_kappa3<-Initial_kappa4<-Initial_kappa5<-
  Initial_kappa6<-Initial_kappa7<-Initial_kappa8<-Initial_kappa9<-Initial_kappa10<-
  vector(,M)

for(m in 1:M){
  Initial_kappa1[m]<-rnorm(1,Initial_mu_1[m],Initial_tau_1[m])
  Initial_kappa2[m]<-rnorm(1,Initial_mu_2[m],Initial_tau_2[m])
  Initial_kappa3[m]<-rnorm(1,Initial_mu_3[m],Initial_tau_3[m])
  Initial_kappa4[m]<-rnorm(1,Initial_mu_4[m],Initial_tau_4[m])
  Initial_kappa5[m]<-rnorm(1,Initial_mu_5[m],Initial_tau_5[m])
  Initial_kappa6[m]<-rnorm(1,Initial_mu_6[m],Initial_tau_6[m])
  Initial_kappa7[m]<-rnorm(1,Initial_mu_7[m],Initial_tau_7[m])
  Initial_kappa8[m]<-rnorm(1,Initial_mu_8[m],Initial_tau_8[m])
  Initial_kappa9[m]<-rnorm(1,Initial_mu_9[m],Initial_tau_9[m])
  Initial_kappa10[m]<-rnorm(1,Initial_mu_10[m],Initial_tau_10[m])
}

#Sample group parameters from the informed prior
Final_mu_1<-rnorm(M,-0.57,0.17*4)
Final_mu_2<-rnorm(M,1.35,0.09*4)
Final_mu_3<-rnorm(M,3.09,0.06*4)
Final_mu_4<-rnorm(M,2.72,0.03*4)
Final_mu_5<-rnorm(M,0.71,0.14*4)
Final_mu_6<-rnorm(M,1.27,0.15*4)
Final_mu_7<-rnorm(M,1.51,0.11*4)
Final_mu_8<-rnorm(M,2.59,0.17*4)
Final_mu_9<-rnorm(M,-3.25,0.39)
Final_mu_10<-rnorm(M,-4.02,0.36)

Final_tau_1<-abs(rnorm(M,1.02,0.15*4))
Final_tau_2<-abs(rnorm(M,0.53,0.07*4))
Final_tau_3<-abs(rnorm(M,0.41,0.05*4))
Final_tau_4<-abs(rnorm(M,0.16,0.04*4))
Final_tau_5<-abs(rnorm(M,0.39,0.19*4))
Final_tau_6<-abs(rnorm(M,0.73,0.14*4))
Final_tau_7<-abs(rnorm(M,0.33,0.14*4))
Final_tau_8<-abs(rnorm(M,0.97,0.15*4))
Final_tau_9<-abs(rnorm(M,0.86,0.36))
Final_tau_10<-abs(rnorm(M,0.39,0.28))

Final_kappa1<-Final_kappa2<-Final_kappa3<-Final_kappa4<-Final_kappa5<-
  Final_kappa6<-Final_kappa7<-Final_kappa8<-Final_kappa9<-Final_kappa10<-
  vector(,M)

for(m in 1:M){
  Final_kappa1[m]<-rnorm(1,Final_mu_1[m],Final_tau_1[m])
  Final_kappa2[m]<-rnorm(1,Final_mu_2[m],Final_tau_2[m])
  Final_kappa3[m]<-rnorm(1,Final_mu_3[m],Final_tau_3[m])
  Final_kappa4[m]<-rnorm(1,Final_mu_4[m],Final_tau_4[m])
  Final_kappa5[m]<-rnorm(1,Final_mu_5[m],Final_tau_5[m])
  Final_kappa6[m]<-rnorm(1,Final_mu_6[m],Final_tau_6[m])
  Final_kappa7[m]<-rnorm(1,Final_mu_7[m],Final_tau_7[m])
  Final_kappa8[m]<-rnorm(1,Final_mu_8[m],Final_tau_8[m])
  Final_kappa9[m]<-rnorm(1,Final_mu_9[m],Final_tau_9[m])
  Final_kappa10[m]<-rnorm(1,Final_mu_10[m],Final_tau_10[m])
}

#Combine and reformat draws
g<-expand.grid(x=seq(0,35,.2),idx=1:M)
kappas<-tibble(
  idx=1:M,
  Initial_kappa1=Initial_kappa1,
  Initial_kappa2=Initial_kappa2,
  Initial_kappa3=Initial_kappa3,
  Initial_kappa4=Initial_kappa4,
  Initial_kappa5=Initial_kappa5,
  Initial_kappa6=Initial_kappa6,
  Initial_kappa7=Initial_kappa7,
  Initial_kappa8=Initial_kappa8,
  Initial_kappa9=Initial_kappa9,
  Initial_kappa10=Initial_kappa10,
  Final_kappa1=Final_kappa1,
  Final_kappa2=Final_kappa2,
  Final_kappa3=Final_kappa3,
  Final_kappa4=Final_kappa4,
  Final_kappa5=Final_kappa5,
  Final_kappa6=Final_kappa6,
  Final_kappa7=Final_kappa7,
  Final_kappa8=Final_kappa8,
  Final_kappa9=Final_kappa9,
  Final_kappa10=Final_kappa10
)

#Estimate and summarize PFs
PFS<-
  full_join(g,kappas)%>%
  pivot_longer(
    cols = !c(x,idx),
    names_to = c('stage','kappa'),
    names_pattern = "(.*)_(.*)"
  ) %>% 
  pivot_wider(
    names_from = kappa
  ) %>% 
  mutate(
    CD=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*(1-2^(-(x/exp(kappa1))^exp(kappa5))),    
    WD=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*(1-2^(-(x/exp(kappa2))^exp(kappa6))),       
    CP=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*(1-2^(-(x/exp(kappa3))^exp(kappa7))),       
    HP=inv_logit(kappa9)/2+(1-inv_logit(kappa9)/2-inv_logit(kappa10)/2)*(1-2^(-(x/exp(kappa4))^exp(kappa8))),    
  ) %>% 
  pivot_longer(
    cols = c(CD,CP,WD,HP),
  ) %>% 
  group_by(name,x,stage) %>% 
  summarise(
    m=median(value),
    lb_60=quantile(value,0.20),
    ub_60=quantile(value,0.80),    
    lb_80=quantile(value,0.10),
    ub_80=quantile(value,0.90),    
    lb_90=quantile(value,0.05),
    ub_90=quantile(value,0.95),
    lb_95=quantile(value,0.025),
    ub_95=quantile(value,0.975)
  ) %>% 
  pivot_longer(
    cols =!c(m,x,name,stage),
    names_to = c('bound','CI'),
    names_pattern = "(.*)_(.*)"
  ) %>% 
  pivot_wider(
    names_from = bound,
    values_from = value
  ) %>% 
  filter((name=='CD'&x<15)|(name=='WD'&x<15)|name=='CP'|(name=='HP'&x<20))

#Make and save the plot
PFS$name<-factor(PFS$name,c('CD','WD','CP','HP'))
PFS$stage<-factor(PFS$stage,c('Initial','Final'))
PFS$x[PFS$name=='CD']<-30-PFS$x[PFS$name=='CD']
PFS$x[PFS$name=='WD']<-30+PFS$x[PFS$name=='WD']
PFS$x[PFS$name=='CP']<-30-PFS$x[PFS$name=='CP']
PFS$x[PFS$name=='HP']<-30+PFS$x[PFS$name=='HP']

PFS %>% 
  ggplot()+
  geom_ribbon(aes(x=x,ymin=lb,ymax=ub,alpha=CI,fill=name))+
  geom_line(aes(x=x,y=m,color=name))+
  facet_grid(stage~name,scales = 'free')+
  labs(
    title='Quick - constrained threshold',
    subtitle='Prior predictive distribution',
    x='Temperature (°C)',
    y='P(response="yes")',
    fill='Modality',
    color='Modality',
    alpha='Quantiles',
  )+
  scale_alpha_manual(labels=c('60%','80%','90%','95%'),values=c(.4,.3,.2,.1))+
  scale_color_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  scale_fill_manual(
    values=c('#56B4E9','#E69F00','#0072B2','#D55E00'),
    labels=c('CD','WD','CP','HP')
  )+
  theme_classic() +
  theme(strip.text.x = element_blank())

ggsave('QCpriors.png',units = 'cm',width = 16,height = 10)
# Script to generate figure 7
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(tidyverse)
library(cmdstanr)
library(ggdist)

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
fit_an<-readRDS("results/fits/Gaussian_constrained_threshold_ageing_and_neuropathy_NRS_NI.rds")

#### Correlation ####
#Extract standardized residuals from ageing and neuropathy fit
draws_df <- fit_an$draws(variables = "partial_residuals", format = "draws_df")
draws_long <-draws_df %>%
  select(starts_with("partial_residuals[")) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(
    parameter = as.integer(str_extract(param, "(?<=\\[)\\d+")),
    participant = as.integer(str_extract(param, "(?<=,)\\d+"))
  ) %>% 
  group_by(parameter,participant) %>% 
  summarize(
    m=mean(value),
    lb=quantile(value,.025),
    ub=quantile(value,.975)
    )

df_pairs <- draws_long %>%
  distinct(parameter) %>%
  rename(param_row = parameter) %>%
  crossing(
    draws_long %>% distinct(parameter) %>% rename(param_col = parameter)
  )
draws_long_col<-
  draws_long %>% 
  rename(
    param_col=parameter,
    m_col=m,
    lb_col=lb,
    ub_col=ub,
  )
draws_long_row<-
  draws_long %>% 
  rename(
    param_row=parameter,
    m_row=m,
    lb_row=lb,
    ub_row=ub,
  )
df_plot<-df_pairs %>% 
  full_join(draws_long_row) %>% 
  full_join(draws_long_col) %>% 
  filter(param_row<9,param_col<9,param_row<param_col)
df_plot$param_col<-factor(df_plot$param_col,1:8,labels=c('CDT','WDT','CPT','HPT','CDS','WDS','CPS','HPS'))
df_plot$param_row<-factor(df_plot$param_row,1:8,labels=c('CDT','WDT','CPT','HPT','CDS','WDS','CPS','HPS'))

#Extract correlation coefficients from ageing and neuropathy fit
draws_df <- fit_an$draws(variables = "omega", format = "draws_df")
correlation_draws_long <- draws_df %>%
  select(starts_with("omega[")) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(
    param_row = as.integer(str_extract(param, "(?<=\\[)\\d+")),
    param_col = as.integer(str_extract(param, "(?<=,)\\d+"))
  )

correlation_draws_filtered <- correlation_draws_long %>%
  filter(param_row > param_col,param_row<9,param_col<9)

correlation_draws_filtered$param_row<-factor(correlation_draws_filtered$param_row,labels=c('WDT','CPT','HPT','CDS','WDS','CPS','HPS'))
correlation_draws_filtered$param_col<-factor(correlation_draws_filtered$param_col,labels=c('CDT','WDT','CPT','HPT','CDS','WDS','CPS'))

#Compute posterior mean and 95% CI for coefficients
rhos<-
  correlation_draws_filtered %>% 
  group_by(param_row,param_col) %>% 
  summarise(
    m=sprintf("%.2f",median(value)),
    lb=sprintf("%.2f",quantile(value,0.025)),
    ub=sprintf("%.2f",quantile(value,0.975))
  )

# Build the plot
facet_text <-
  rhos %>%
  mutate(
    label = paste0('ρ = ',m, "\n[", lb, ", ", ub, "]")
  )

corr_plot<-
  df_plot %>% 
  ggplot() +
  geom_pointrange(aes(x=m_col,y=m_row,ymin=lb_row,ymax=ub_row),alpha=.1,size=.1)+
  geom_linerange(aes(x=m_col,y=m_row,xmin=lb_col,xmax=ub_col),alpha=.1)+
  geom_text(
    data = facet_text,
    aes(
      x = 0,
      y = 0,
      label = label
    ),
    inherit.aes = FALSE,
    size = 3
  ) +
  facet_grid(param_row~param_col,scales = 'free') +
  theme_classic() +
  labs(
    title ='Partial correlations between parameters',
    x = "", 
    y = ""
  )+
  xlim(c(-4,4))+
  ylim(c(-4,4))
    
ggsave(
  'plots/corr_plot_mod.jpeg',
  plot = corr_plot,
  width = 18,
  height = 18,
  dpi = 300,
  units = 'cm'
)

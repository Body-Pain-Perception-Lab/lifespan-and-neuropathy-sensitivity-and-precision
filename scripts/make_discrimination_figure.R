# Script to test discrimination accuracy and make discrimination figure
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

#### Set-up environment ####
library(cmdstanr)
library(tidyverse)
library(pROC)

rm(list=ls())

directory=getwd()

#### Load data ####
fit<-readRDS(paste0(directory,"/results/fits/Gaussian_single_subjects.rds"))
d <- read_csv("data/aggregated_results.csv") %>% 
  filter(
    recording_deviates_from_target==0,
    recording_deviates_from_mean==0,
    !is.na(response),
    !subject%in%c(3019) #exclusion based on saturation check plots
  ) %>% 
  arrange(subject)
participant<-unique(d$subject)

participant_info<-read_csv("data/demographics.csv") %>% 
  filter(ID%in%participant) %>% 
  arrange(ID)

#### Extract posterior means ####
posterior_means<-
  fit$summary(c('alpha','beta')) %>% 
  select(variable,mean) %>% 
  extract(
    variable,
    into = c("parameter", "condition", "participant"),
    regex = "([A-Za-z]+)\\[(\\d+),\\s*(\\d+)\\]"
    ) %>% 
  mutate(col_name = paste0(parameter, "_", condition)) %>%
  select(participant,col_name,mean) %>% 
  pivot_wider(
    names_from = col_name,
    values_from = mean
  ) %>% 
  rename(
    cdt=alpha_1,
    wdt=alpha_2,
    cpt=alpha_3,
    hpt=alpha_4,
    l_cds=beta_1,
    l_wds=beta_2,
    l_cps=beta_3,
    l_hps=beta_4
  )
posterior_means$age<-participant_info$age
posterior_means$status<-participant_info$ID>=5000

#### Construct models and assess their discrimination performance ####
model_list <- list(
  CDT   = status ~ cdt + cdt:age,
  CDS   = status ~ l_cds + l_cds:age,
  `CDT&CDS`    = status ~ cdt + l_cds + cdt:age + l_cds:age,
  
  WDT   = status ~ wdt + wdt:age,
  WDS   = status ~ l_wds + l_wds:age,
  `WDT&WDS`    = status ~ wdt + l_wds + wdt:age + l_wds:age,
  
  CPT   = status ~ cpt + cpt:age,
  CPS   = status ~ l_cps + l_cps:age,
  `CPT&CPS`    = status ~ cpt + l_cps + cpt:age + l_cps:age,
  
  HPT   = status ~ hpt + hpt:age,
  HPS   = status ~ l_hps + l_hps:age,
  `HPT&HPS` = status ~ hpt + l_hps + hpt:age + l_hps:age,
  
  ALL = status ~ 
    cdt + wdt + cpt + hpt +
    l_cds + l_wds + l_cps + l_hps +
    cdt:age + wdt:age + cpt:age + hpt:age +
    l_cds:age + l_wds:age + l_cps:age + l_hps:age
)

results <- data.frame(
  model = character(),
  auc = numeric(),
  ci_low = numeric(),
  ci_high = numeric(),
  stringsAsFactors = FALSE
)
roc_store<-list()

for (m in names(model_list)) {
  
  fit <- glm(
    model_list[[m]],
    data = posterior_means,
    family = binomial
  )
  
  prob <- predict(fit, type = "response")
  
  roc_obj <- roc(
    response = posterior_means$status,
    predictor = prob
  )
  
  roc_store[[m]] <- roc_obj
  
  set.seed(1)
  ci_vals <- ci.auc(
    roc_obj,
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    conf.level = 0.95
  )
  
  results <- rbind(
    results,
    data.frame(
      model = m,
      auc = as.numeric(auc(roc_obj)),
      ci_low = ci_vals[1],
      ci_high = ci_vals[3],
      conf.level = 0.95
      )
  )
}

plot_data <- 
  results %>%
  mutate(group =
           case_when(
             grepl("^CDT|^CDS", model) ~ "CD",
             grepl("^WDT|^WDS", model) ~ "WD",
             grepl("^CPT|^CPS", model) ~ "CP",
             grepl("^HPT|^HPS", model) ~ "HP",
             grepl("^ALL", model) ~ "ALL"
           )
  )
plot_data$group<-factor(plot_data$group,c('ALL','CD','WD','CP','HP'))
discrimination_plot<-
  plot_data %>% 
  ggplot()+
  geom_pointrange(aes(x=model,y=auc,ymin=ci_low,ymax=ci_high,color=group))+
  geom_hline(aes(yintercept = .5),linetype='dashed')+
  geom_hline(aes(yintercept = 1),linetype='dashed')+
  scale_color_manual(
    values = c('#6C3BAA','#56B4E9', '#E69F00', '#0072B2', '#D55E00'),
    labels = c('All','CD', 'WD', 'CP', 'HP')
  ) +  
  labs(
    color='',
    x='Model',
    y='AUROCC'
    )+
  theme_classic()+
  theme(legend.position = 'bottom')

ggsave(
  'plots/discrimination.jpeg',
  plot = discrimination_plot,
  width = 18,
  height = 10,
  dpi = 300,
  units = 'cm'
)

set.seed(1)
cdt_cd<-
  roc.test(
  roc_store[[1]],
  roc_store[[3]],
  method = "bootstrap",
  boot.n = 5000,
  boot.stratified = TRUE,
  alternative='less'
  )
set.seed(1)
cds_cd<-
  roc.test(
    roc_store[[2]],
    roc_store[[3]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )

set.seed(1)
wdt_wd<-
  roc.test(
    roc_store[[4]],
    roc_store[[6]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )
set.seed(1)
wds_wd<-
  roc.test(
    roc_store[[5]],
    roc_store[[6]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )

set.seed(1)
cpt_cp<-
  roc.test(
    roc_store[[7]],
    roc_store[[9]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )
set.seed(1)
cps_cp<-
  roc.test(
    roc_store[[8]],
    roc_store[[9]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )

set.seed(1)
hpt_hp<-
  roc.test(
    roc_store[[10]],
    roc_store[[12]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )
set.seed(1)
hps_hp<-
  roc.test(
    roc_store[[11]],
    roc_store[[12]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )

set.seed(1)
all_cd<-
  roc.test(
    roc_store[[3]],
    roc_store[[13]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )
set.seed(1)
all_wd<-
  roc.test(
    roc_store[[6]],
    roc_store[[13]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )
set.seed(1)
all_cp<-
  roc.test(
    roc_store[[9]],
    roc_store[[13]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )
set.seed(1)
all_hp<-
  roc.test(
    roc_store[[12]],
    roc_store[[13]],
    method = "bootstrap",
    boot.n = 5000,
    boot.stratified = TRUE,
    alternative='less'
  )

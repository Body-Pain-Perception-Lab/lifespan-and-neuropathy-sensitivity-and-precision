# Script to summarize the sample demographics data
# Author: Arthur S. Courtin  
# License: MIT (see LICENSE file) 

demographics %>% 
  filter(!ID%in%participant_to_exclude$ID,ID!=3019) %>% 
  mutate(status=ID>=5000) %>% 
  group_by(status) %>% 
  summarise(
    n=sum(!is.na(gender)),
    n_M=sum(gender=='M'),
    n_F=sum(gender=='F'),
    
    min_age=min(age),
    max_age=max(age)
  )

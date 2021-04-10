### COnsumption by size 

##df_p1 <- read.csv("data/Unused Data/size_experiment/raw/Purple_size_data.csv")

df_p1 <- df_p1 %>% 
  mutate(urchin = "p")

str(df_p1)
head(df_p1)

df_p2 <- read.csv("data/Unused Data/size_experiment/raw/Purple_consuption_data.csv")

df_p3 <-df_p2 %>%
  filter(urchin_size_cat != 'control') %>% 
  mutate(urchin = "p")%>% 
  left_join(df_p1) %>% 
  select(-date) %>% 
  select (-days_starved) %>% 
  select (-time_in) %>% 
  select (-time_out) %>% 
  #select (-kelp_per_urchin) %>%
  #select (-kelp_consumed) %>% 
  mutate(level= "B") %>% 
  mutate (time_ran = 48) %>% 
  mutate( round= week_no) %>% 
  select(-week_no) %>% 
  mutate(size_class= urchin_size_cat) %>% 
  select (-urchin_size_cat) %>% 
  mutate (urchin_density =urchin_dens) %>% 
  select (-urchin_dens) %>% 
  mutate(kelp_consumed=(kelp_in-kelp_out)/48) %>% 
  group_by(trial_id) %>% 
  mutate(per_cap=kelp_consumed/7)



plot (df_p3$size, df_p3$per_cap) 


lmsize = lm(per_cap~size, data = df_p3) #Create the linear regression
summary(lmsize) #Review the results

#consumption=a+size*b

b<-.00089
a<-.003482
size_ln<-(df_p3$size)

p4<-df_p3 %>% 
  mutate(regression=a+b*size)

plot (df_p3$size, df_p3$per_cap, col="green") 
lines(p4$size, p4$regression, col="blue" )

# need to fit other ones to it and do AIC to prove that linear regression was the best choice

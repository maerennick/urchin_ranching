
#### Trying a new appraoch to the model 

gamma_1<-0.99 #100% of 1cm urchins will transition into the 2cm size class within a year
gamma_2<-0.8
gamma_3<-0.7
gamma_4<-0.6
gamma_5<-0.5
gamma_6<-0.4

#Death -> add a death function here as well 

Z_1<-.9
Z_2<-.7
Z_3<-.5
Z_4<-.3
Z_5<-.1
Z_6<-.1


#other parameters

Kmax<-50000
R<-0.5
u<-0
r<- (.36*0.5*10) #survivability*% of the reproducing population presumed to be female #(recruitment function)



b<-.00089
a<-.003482

c1<- (1.5*24*365*a+b)
c2<- (2.5*24*365*a+b)
c3<- (3.5*24*365*a+b)
c4<- (4.5*24*365*a+b)
c5<- (6.5*24*365*a+b)
c6<- (7.5*24*365*a+b)

tset <- seq(from = 0, to = 100, length.out = 100)





k.simu <- NaN*tset; k.simu[1] <- 10000
u1.simu <- NaN*tset; u1.simu[1] <- 100
u2.simu <- NaN*tset; u2.simu[1] <- 200
u3.simu <- NaN*tset; u3.simu[1] <- 500
u4.simu <- NaN*tset; u4.simu[1] <- 500
u5.simu <- NaN*tset; u5.simu[1] <- 100
u6.simu <- NaN*tset; u6.simu[1] <- 50


for(i in 2:100){
  # calculate change in time
  dt <- tset[i] - tset[i-1]
  
  # store dummy variables
  k <- k.simu[i-1]
  u1 <- u1.simu[i-1]
  u2 <- u2.simu[i-1]
  u3 <- u3.simu[i-1]
  u4 <- u4.simu[i-1]
  u5 <- u5.simu[i-1]
  u6 <- u6.simu[i-1]
  
  
  # calculate change in population size
  
  
  dk <- ((1-R*(1-(k/Kmax)))*k)*dt
  du1 <- (((1-gamma_1*exp(-Z_1)*u1))+ (r*(u4+u5+u6)))*dt
  du2 <-  ((gamma_1*exp(-Z_1)*u1)+(((1-gamma_2)*exp(-Z_2)*u2)))*dt
  du3 <-  ((gamma_2*exp(-Z_2)*u2)+(((1-gamma_3)*exp(-Z_3)*u3)))*dt
  du4 <-  ((gamma_3*exp(-Z_3)*u3)+(((1-gamma_4)*exp(-Z_4)*u4)))*dt
  du5 <-  ((gamma_4*exp(-Z_4)*u4)+(((1-gamma_5)*exp(-Z_5)*u5)))*dt
  du6 <-  ((gamma_5*exp(-Z_5)*u5)+(((1-gamma_6)*exp(-Z_6)*u6)))*dt
  
  
  # calculate total population size
  u1.simu[i]<-du1
  u2.simu[i]<-du2
  u3.simu[i]<-du3
  u4.simu[i]<-du4
  u5.simu[i]<-du5
  u6.simu[i]<-du6
  k.simu [i]<-k.simu[i-1]+ dk

}


k.simu

#u1.simu
u2.simu
#u3.simu
#u4.simu
#u5.simu
#u6.simu

plot(x = tset, y = u1.simu, col = 'black', type='p', lwd = 2, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,30),
     ylim=c(0,2000))
lines(x = tset, y = u2.simu, col = 'blue', lwd = 2)
lines(x = tset, y = u3.simu, col = 'green', lwd = 2,lty = 2)
lines(x = tset, y = u4.simu, col = 'dark blue', lwd = 2,lty = 2)
lines(x = tset, y = u5.simu, col = 'light blue', lwd = 2,lty = 2)
lines(x = tset, y = u6.simu, col = 'dark green', lwd = 2,lty = 2)
legend(x = max(30)*0.849, y = max(2000), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('black', 'blue', 'green','dark blue', 'light blue', 'dark green' ))


plot(x = tset, y = k.simu, col = 'black', type='p', lwd = 2, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,30),
     ylim=c(0,1000000))




#### Adding in real world dynamics (bring in kelp dynamics and add recruitment)#####################################################################################################################################################################################################

## Von Bertalanffy Equation 
#Test diameter as a product of age

#Din<- 63.38
#L<- .327 #growth rate constant
#a0<-0

diameter <- function(a) {
  D <- Din*(1-exp(-L*(a-a0)))
  return(D)
}

age<-c(0,0.25,0.5,0.75,1, 1.25,1.5,1.75,2, 2.25,2.5,2.75,3, 3.25,3.5,3.75,4, 4.25,4.5,4.75,5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7)

diameter(age)

plot(diameter)


gamma_1<-.99 #100% of 1cm urchins will transition into the 2cm size class within a year
gamma_2<-.9
gamma_3<-0.7
gamma_4<-0.5
gamma_5<-0.3
gamma_6<-0.1

#Death -> add a death function here as well 

library(tidyverse)
library(ggpubr)

# Get size data

urc <- read.csv("data/LTE_Urchin_All_Years_20190611.csv", header = T) %>% # LTER dataset: urchin size and frequency data collected from 5 sites between 2008-2019 in the Santa Barbara Channel. #They only take this data once at each site? 
  filter(TREATMENT == "CONTROL") %>% select(YEAR, MONTH, DATE, SITE, TRANSECT, SIZE, COUNT, COMMON_NAME) %>% 
  rename_all(tolower) %>% 
  group_by(year, month, date, site, transect, common_name, size) %>% 
  summarize(count = sum(count)) %>% #What does this change other than moving it? 
  group_by(year, month, date, site, transect, common_name, size) %>%
  complete(count = full_seq(1:count,1))%>% #transferring each count to one, and seperating each urchin into its own line in the data. 
  select(-count) %>% #removing the count column
  ungroup() %>%
  mutate(size = as.numeric(size))%>% 
  filter(common_name== "Purple Urchin") %>%
  filter(site=="NAPL")

urc_count<-urc %>%
  group_by(year, size) %>% 
  summarise(n = n()) %>%
  ungroup() %>% 
  group_by(size) %>% 
  mutate(avg=(mean(n)))
  


Z_1<-.9
Z_2<-.5
Z_3<-.3
Z_4<-.3
Z_5<-.1
Z_6<-.1




#other parameters

K<-1500
R<-.5
u<-0
r<- 1000 #survivability*% of the reproducing population presumed to be female #(recruitment function) # 1000 mew babies over 60 m^2 (aka transect coverage)
K_u<-600
Ro<-6
a_1<-1.25
a_2<-0.05
b_1<-0.15
b_2<-0.05
alpha<-1


b<-.0032674
a<-.009645

c1<- (1.5*24*365*a)*0.001
c2<- (2.5*24*365*a)*0.001
c3<- (3.5*24*365*a)*0.001
c4<- (4.5*24*365*a)*0.001
c5<- (6.5*24*365*a)*0.001
c6<- (7.5*24*365*a)*0.001

tset <- seq(from = 0, to = 100, length.out = 100)


dr.simu <- NaN*tset; dr.simu[1] <- 1.6
dd.simu <- NaN*tset; dr.simu[1] <- 0.2
du.simu <- NaN*tset; dr.simu[1] <- 0.2


k.simu <- NaN*tset; k.simu[1] <- 1000
u1.simu <- NaN*tset; u1.simu[1] <- 20
u2.simu <- NaN*tset; u2.simu[1] <- 15
u3.simu <- NaN*tset; u3.simu[1] <- 30
u4.simu <- NaN*tset; u4.simu[1] <- 15
u5.simu <- NaN*tset; u5.simu[1] <- 10
u6.simu <- NaN*tset; u6.simu[1] <- 5


 
for(i in 2:100){
  # calculate change in time
  dt <- tset[i] - tset[i-1]
  
  # store dummy variables
  
  dr<-dr.simu[i-1]
  dd<-dd.simu[i-1]
  du<-du.simu[i-1]
    
  k <- k.simu[i-1]
  u1 <- u1.simu[i-1]
  u2 <- u2.simu[i-1]
  u3 <- u3.simu[i-1]
  u4 <- u4.simu[i-1]
  u5 <- u5.simu[i-1]
  u6 <- u6.simu[i-1]
  
  
  # calculate change in population size
  
  
  #(R*k/((1+((R-1)/(K-(u1*c1+u2*c2+u3*c3+u4*c4+u5*c5+u6*c6))))*k)*dt
  #dk<-((k-((u1*c1+u2*c2+u3*c3+u4*c4+u5*c5+u6*c6)*k))*((exp(c))/(1+d*(k-((u1*c1+u2*c2+u3*c3+u4*c4+u5*c5+u6*c6)*k)))))*dt
  
  
  dr<-a_1+(a_2/2)*(1-cos(pi*2*tset[i]))*dt
  dd<-b_1-b_2*sin(pi*2*tset[i])*dt
  du<-(c1+0.0032674)*u1+
    (c2+0.0032674)*u2+
    (c3+0.0032674)*u3+
    (c4+0.0032674)*u4+
    (c5+0.0032674)*u5+
    (c6+0.0032674)*u6
    
    
    
  dk<- (((1-(k/K))*(dr-dd))*k)-du
  
  
  du1 <- ((Ro*(u4+u5+u6)*(1-((u4+u5+u6)/K_u)))-(gamma_1*exp(-Z_1)*u1))*dt
  du2 <-  ((gamma_1*exp(-Z_1)*u1)+(((1-gamma_2)*exp(-Z_2)*u2)))*dt
  du3 <-  ((gamma_2*exp(-Z_2)*u2)+(((1-gamma_3)*exp(-Z_3)*u3)))*dt
  du4 <-  ((gamma_3*exp(-Z_3)*u3)+(((1-gamma_4)*exp(-Z_4)*u4)))*dt
  du5 <-  ((gamma_4*exp(-Z_4)*u4)+(((1-gamma_5)*exp(-Z_5)*u5)))*dt
  du6 <-  ((gamma_5*exp(-Z_5)*u5)+(((1-gamma_6)*exp(-Z_6)*u6)))*dt
  
  
  # calculate total population size
  u1.simu[i]<-du1
  u2.simu[i]<-du2
  u3.simu[i]<-du3
  u4.simu[i]<-du4
  u5.simu[i]<-du5
  u6.simu[i]<-du6
  k.simu [i]<-k.simu[i-1]+ dk
  du.simu [i]<-du
  dr.simu [i]<-dr
  dd.simu [i]<-dd
  
}


u1*c1+u2*c2+u3*c3+u4*c4+u5*c5+u6*c6

k.simu

dr.simu
dd.simu

u1.simu
u2.simu
#u3.simu
#u4.simu
#u5.simu
#u6.simu

plot(x = tset, y = u1.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,10),
     ylim=c(0,600))
lines(x = tset, y = u2.simu, col = 'blue', lwd = 2)
lines(x = tset, y = u3.simu, col = 'green', lwd = 2,lty = 2)
lines(x = tset, y = u4.simu, col = 'dark blue', lwd = 2,lty = 2)
lines(x = tset, y = u5.simu, col = 'light blue', lwd = 2,lty = 2)
lines(x = tset, y = u6.simu, col = 'dark green', lwd = 2,lty = 2)
legend(x = max(30)*0.849, y = max(2000), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('black', 'blue', 'green','dark blue', 'light blue', 'dark green' ))


plot(x = tset, y = k.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,2000))


################################################
######### Add in fishing pressure


Z_1<-.9
Z_2<-.5
Z_3<-.3
Z_4<-(.3+0.4)
Z_5<-(.1+0.6) #(death rate + fishing pressure)
Z_6<-(.1+0.6) #(death rate + fishing pressure)



for(i in 2:100){
  # calculate change in time
  dt <- tset[i] - tset[i-1]
  
  # store dummy variables
  
  dr<-dr.simu[i-1]
  dd<-dd.simu[i-1]
  du<-du.simu[i-1]
  
  k <- k.simu[i-1]
  u1 <- u1.simu[i-1]
  u2 <- u2.simu[i-1]
  u3 <- u3.simu[i-1]
  u4 <- u4.simu[i-1]
  u5 <- u5.simu[i-1]
  u6 <- u6.simu[i-1]
  
  
  dr<-a_1+(a_2/2)*(1-cos(pi*2*tset[i]))*dt
  dd<-b_1-b_2*sin(pi*2*tset[i])*dt
  du<-(c1+0.0032674)*u1+
    (c2+0.0032674)*u2+
    (c3+0.0032674)*u3+
    (c4+0.0032674)*u4+
    (c5+0.0032674)*u5+
    (c6+0.0032674)*u6
  
  dk<- (((1-(k/K))*(dr-dd))*k)-du
  
  
  du1 <- ((Ro*(u4+u5+u6)*(1-((u4+u5+u6)/K_u)))-(gamma_1*exp(-Z_1)*u1))*dt
  du2 <-  ((gamma_1*exp(-Z_1)*u1)+(((1-gamma_2)*exp(-Z_2)*u2)))*dt
  du3 <-  ((gamma_2*exp(-Z_2)*u2)+(((1-gamma_3)*exp(-Z_3)*u3)))*dt
  du4 <-  ((gamma_3*exp(-Z_3)*u3)+(((1-gamma_4)*exp(-Z_4)*u4)))*dt
  du5 <-  ((gamma_4*exp(-Z_4)*u4)+(((1-gamma_5)*exp(-Z_5)*u5)))*dt
  du6 <-  ((gamma_5*exp(-Z_5)*u5)+(((1-gamma_6)*exp(-Z_6)*u6)))*dt
  
  
  # calculate total population size
  u1.simu[i]<-du1
  u2.simu[i]<-du2
  u3.simu[i]<-du3
  u4.simu[i]<-du4
  u5.simu[i]<-du5
  u6.simu[i]<-du6
  k.simu [i]<-k.simu[i-1]+ dk
  du.simu [i]<-du
  dr.simu [i]<-dr
  dd.simu [i]<-dd
  
}

plot(x = tset, y = u1.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,10),
     ylim=c(0,600))
lines(x = tset, y = u2.simu, col = 'blue', lwd = 2)
lines(x = tset, y = u3.simu, col = 'green', lwd = 2,lty = 2)
lines(x = tset, y = u4.simu, col = 'dark blue', lwd = 2,lty = 2)
lines(x = tset, y = u5.simu, col = 'light blue', lwd = 2,lty = 2)
lines(x = tset, y = u6.simu, col = 'dark green', lwd = 2,lty = 2)
legend(x = max(30)*0.849, y = max(2000), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('black', 'blue', 'green','dark blue', 'light blue', 'dark green' ))


plot(x = tset, y = k.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,2000))




#####################sensitivity 







par_names<- c("K |change = abschg", "Ro |change = pctchg", "K_u |change = abschg", "gamma_1 |change = pctchg", "Z_1|change = pctchg", "gamma_2|change = pctchg", "Z_2|change = pctchg", "gamma_3|change = pctchg", "Z_3|change = pctchg", "gamma_4|change = pctchg", "Z_4|change = pctchg", "gamma_5|change = pctchg", "Z_5|change = pctchg", "gamma_6|change = pctchg", "Z_6|change = pctchg")
  
  
  "cn2.hru|change = pctchg" = - 5)


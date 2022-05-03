#### Mae Rennick
##EEMB595CA
###FINAL Project Code
### Integrative approach to quantifying the conservation potential of urchin removal for kelp restoration



###############Parameters 

K<-400 #Max kelp 
u<-0 #harvest of kelp 
K_u<-400 #max urchins
Ro<-4.5 #growth rate 
a_1<-1.25# unknown regression parameter
a_2<-0.05 #unknown regression parameter
b_1<-.15 #unknown regression parameter
b_2<-0.05 #unknown regression parameter


gamma_1<-.99 #probability of surviving urchins in the N1 size class transitioning to the N2 size class
gamma_2<-.9 #probability of surviving urchins in the N2 size class transitioning to the N3 size class
gamma_3<-0.7 #probability of surviving urchins in the N3 size class transitioning to the N4 size class
gamma_4<-0.5 #probability of surviving urchins in the N4 size class transitioning to the N5 size class
gamma_5<-0.3 #probability of surviving urchins in the N5 size class transitioning to the N6 size class
gamma_6<-0.1 #probability of surviving urchins in the N6 size class dying


b<-.0032674 #consumption parameter from Rennick et al. in review (y intercept)
a<-.009645 #consumption parameter from Rennick et al. in review (slope)

c1<- (1.5*24*365*a)*0.001  #Consumption of urchins in the N1 size class kg/year
c2<- (2.5*24*365*a)*0.001  #Consumption of urchins in the N2 size class kg/year
c3<- (3.5*24*365*a)*0.001  #Consumption of urchins in the N3 size class kg/year
c4<- (4.5*24*365*a)*0.001  #Consumption of urchins in the N4 size class kg/year
c5<- (6.5*24*365*a)*0.001  #Consumption of urchins in the N5 size class kg/year
c6<- (7.5*24*365*a)*0.001  #Consumption of urchins in the N6 size class kg/year

tset <- seq(from = 0, to = 100, length.out = 100) #time steps


dr.simu <- NaN*tset; dr.simu[1] <- 1.6 #empty vector for dr simulation and initial value
dd.simu <- NaN*tset; dr.simu[1] <- 0.2 #empty vector for dd simulation and initial value
du.simu <- NaN*tset; dr.simu[1] <- 0.2 #empty vector for du simulation and initial value for kelp consumption


k.simu <- NaN*tset; k.simu[1] <- 100 #empty vector for k simulation and initial value for kelp
u1.simu <- NaN*tset; u1.simu[1] <- 20 #empty vector for u1 simulation and initial value of N1 urchins
u2.simu <- NaN*tset; u2.simu[1] <- 15 #empty vector for u2 simulation and initial value of N2 urchins
u3.simu <- NaN*tset; u3.simu[1] <- 30 #empty vector for u3 simulation and initial value of N3 urchins
u4.simu <- NaN*tset; u4.simu[1] <- 15 #empty vector for u4 simulation and initial value of N4 urchins
u5.simu <- NaN*tset; u5.simu[1] <- 10 #empty vector for u5 simulation and initial value of N5 urchins
u6.simu <- NaN*tset; u6.simu[1] <- 5 #empty vector for u6 simulation and initial value of N6 urchins


############### MODEL with No fishing #################################################################################################################################################################################################

Z_1<-.9 #annual probability of death+ fishing pressure for N1
Z_2<-.5 #annual probability of death+ fishing pressure for N2
Z_3<-.3 #annual probability of death+ fishing pressure for N3
Z_4<-.3 #annual probability of death+ fishing pressure for N4
Z_5<-.1 #annual probability of death+ fishing pressure for N5
Z_6<-.1 #annual probability of death+ fishing pressure for N6





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
  

  # calculate change in population size for each state variable 

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
  k.simu[i]<-k.simu[i-1]+ dk
  du.simu [i]<-du
  dr.simu [i]<-dr
  dd.simu [i]<-dd
  
}

#Plot kelp biomass through time 

plot(x = tset, y = k.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,300))

#Plot size-specific urchin abundance through time 

plot(x = tset, y = u1.simu, col = '#43225c', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,50),
     ylim=c(0,500), main= "Modeling Size-Structured Urchin Abundance")
lines(x = tset, y = u2.simu, col = '#897295', lwd = 4)
lines(x = tset, y = u3.simu, col = '#8284bc', lwd = 4)
lines(x = tset, y = u4.simu, col = '#cecee9', lwd = 4)
lines(x = tset, y = u5.simu, col = '#695973', lwd = 4)
lines(x = tset, y = u6.simu, col = '#c7b6c8', lwd = 4)
legend(x = max(2)*0.849, y = max(500), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('#43225c', '#897295', '#8284bc','#cecee9', '#695973', '#c7b6c8' ), horiz=TRUE, bty = "n")


library(dplyr)

library(tidyverse)

## MODEL with mild fishing   #################################################################################################################################################################################################

k.simu_2 <- NaN*tset; k.simu_2[1] <- 100


Z_1_2<-.9 #annual probability of death+ fishing pressure for N1
Z_2_2<-.5 #annual probability of death+ fishing pressure for N2
Z_3_2<-.3 #annual probability of death+ fishing pressure for N3
Z_4_2<-(.3+0.4) #(death rate + fishing pressure)
Z_5_2<-(.1+0.4) #(death rate + fishing pressure)
Z_6_2<-(.1+0.4) #(death rate + fishing pressure)



for(i in 2:100){
  # calculate change in time
  dt <- tset[i] - tset[i-1]
  
  # store dummy variables
  
  dr<-dr.simu[i-1]
  dd<-dd.simu[i-1]
  du<-du.simu[i-1]
  
  k <- k.simu_2[i-1]
  u1 <- u1.simu[i-1]
  u2 <- u2.simu[i-1]
  u3 <- u3.simu[i-1]
  u4 <- u4.simu[i-1]
  u5 <- u5.simu[i-1]
  u6 <- u6.simu[i-1]
  
  
  # calculate change in population size
  
  
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
  du2 <-  ((gamma_1*exp(-Z_1_2)*u1)+(((1-gamma_2)*exp(-Z_2_2)*u2)))*dt
  du3 <-  ((gamma_2*exp(-Z_2_2)*u2)+(((1-gamma_3)*exp(-Z_3_2)*u3)))*dt
  du4 <-  ((gamma_3*exp(-Z_3_2)*u3)+(((1-gamma_4)*exp(-Z_4_2)*u4)))*dt
  du5 <-  ((gamma_4*exp(-Z_4_2)*u4)+(((1-gamma_5)*exp(-Z_5_2)*u5)))*dt
  du6 <-  ((gamma_5*exp(-Z_5_2)*u5)+(((1-gamma_6)*exp(-Z_6_2)*u6)))*dt
  
  
  # calculate total population size
  
  u1.simu[i]<-du1
  u2.simu[i]<-du2
  u3.simu[i]<-du3
  u4.simu[i]<-du4
  u5.simu[i]<-du5
  u6.simu[i]<-du6
  k.simu_2 [i]<-k.simu_2[i-1]+ dk
  du.simu [i]<-du
  dr.simu [i]<-dr
  dd.simu [i]<-dd
  
}

#Plot kelp biomass through time 

plot(x = tset, y = k.simu_2, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,2000))


#Plot size-specific urhcin abundance through time

plot(x = tset, y = u1.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,100),
     ylim=c(0,600))
lines(x = tset, y = u2.simu, col = 'blue', lwd = 2)
lines(x = tset, y = u3.simu, col = 'green', lwd = 2,lty = 2)
lines(x = tset, y = u4.simu, col = 'dark blue', lwd = 2,lty = 2)
lines(x = tset, y = u5.simu, col = 'light blue', lwd = 2,lty = 2)
lines(x = tset, y = u6.simu, col = 'dark green', lwd = 2,lty = 2)
legend(x = max(30)*0.849, y = max(2000), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('black', 'blue', 'green','dark blue', 'light blue', 'dark green' ))



## MODEL with heavy fishing   #################################################################################################################################################################################################
k.simu_3 <- NaN*tset; k.simu_3[1] <- 100

Z_1_3<-.9 #annual probability of death+ fishing pressure for N1
Z_2_3<-.5 #annual probability of death+ fishing pressure for N2
Z_3_3<-(.3+0.69) #(death rate + fishing pressure)
Z_4_3<-(.3+0.69) #(death rate + fishing pressure)
Z_5_3<-(.1+0.89) #(death rate + fishing pressure)
Z_6_3<-(.1+0.89) #(death rate + fishing pressure)



for(i in 2:100){
  # calculate change in time
  dt <- tset[i] - tset[i-1]
  
  # store dummy variables
  
  dr<-dr.simu[i-1]
  dd<-dd.simu[i-1]
  du<-du.simu[i-1]
  
  k <- k.simu_3[i-1]
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
  du2 <-  ((gamma_1*exp(-Z_1_3)*u1)+(((1-gamma_2)*exp(-Z_2_3)*u2)))*dt
  du3 <-  ((gamma_2*exp(-Z_2_3)*u2)+(((1-gamma_3)*exp(-Z_3_3)*u3)))*dt
  du4 <-  ((gamma_3*exp(-Z_3_3)*u3)+(((1-gamma_4)*exp(-Z_4_3)*u4)))*dt
  du5 <-  ((gamma_4*exp(-Z_4_3)*u4)+(((1-gamma_5)*exp(-Z_5_3)*u5)))*dt
  du6 <-  ((gamma_5*exp(-Z_5_3)*u5)+(((1-gamma_6)*exp(-Z_6_3)*u6)))*dt
  
  
  # calculate total population size
  u1.simu[i]<-du1
  u2.simu[i]<-du2
  u3.simu[i]<-du3
  u4.simu[i]<-du4
  u5.simu[i]<-du5
  u6.simu[i]<-du6
  k.simu_3 [i]<-k.simu_3[i-1]+ dk
  du.simu [i]<-du
  dr.simu [i]<-dr
  dd.simu [i]<-dd
  
}

#Plot kelp biomass through time 

plot(x = tset, y = k.simu_3, col = 'dark green', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,50),
     ylim=c(0,2000))

#Plot size-specific urhcin abundance through time

plot(x = tset, y = u1.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,40),
     ylim=c(0,200))
lines(x = tset, y = u2.simu, col = 'blue', lwd = 2)
lines(x = tset, y = u3.simu, col = 'green', lwd = 2,lty = 2)
lines(x = tset, y = u4.simu, col = 'dark blue', lwd = 2,lty = 2)
lines(x = tset, y = u5.simu, col = 'light blue', lwd = 2,lty = 2)
lines(x = tset, y = u6.simu, col = 'dark green', lwd = 2,lty = 2)
legend(x = max(30)*0.849, y = max(2000), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('black', 'blue', 'green','dark blue', 'light blue', 'dark green' ))






##################################### Plot with kelp biomass at each level of fishing pressure 



plot(x = tset, y = k.simu, col = '#006a4e', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='kelp abundance (kg)',
     xlim=c(0,100),
     ylim=c(0,600))
lines(x = tset, y = k.simu_2, col = '#5ca08e', lwd = 4)
lines(x = tset, y = k.simu_3, col = '#b8d5cd', lwd = 4)
legend(x = max(100)*0.7, y = max(600), legend = c('No fishing','Fishing', 'High Fishing'),lwd = 4, col =c('#006a4e', '#5ca08e', '#b8d5cd'))



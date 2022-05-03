##### Updated model
## Mae Rennick
### 03-03-22



###############Parameters 

K<-4000 #Max kelp (based on a 100 M^2 kelp forrest? in kg) 
u<-0 #harvest of kelp 
K_u<-400 #max urchins
Ro<-2 #recruitment rate per reproductive individual (urchin)
a_1<-1.25# unknown regression parameter
a_2<-0.05 #unknown regression parameter
b_1<-.15 #unknown regression parameter
b_2<-0.05 #unknown regression parameter


i_gamma_1<-.99 #probability of surviving urchins in the N1 size class transitioning to the N2 size class
i_gamma_2<-.9 #probability of surviving urchins in the N2 size class transitioning to the N3 size class
i_gamma_3<-0.7 #probability of surviving urchins in the N3 size class transitioning to the N4 size class
i_gamma_4<-0.5 #probability of surviving urchins in the N4 size class transitioning to the N5 size class
i_gamma_5<-0.3 #probability of surviving urchins in the N5 size class transitioning to the N6 size class
i_gamma_6<-0.1 #probability of surviving urchins in the N6 size class dying


b<-.0032674 #consumption parameter from Rennick et al. 2022 (y intercept)
a<-.009645 #consumption parameter from Rennick et al. 2022 (slope)

c1<- (1.5*24*365*a)*0.001  #Consumption of urchins in the N1 size class kg/year
c2<- (2.5*24*365*a)*0.001  #Consumption of urchins in the N2 size class kg/year
c3<- (3.5*24*365*a)*0.001  #Consumption of urchins in the N3 size class kg/year
c4<- (4.5*24*365*a)*0.001  #Consumption of urchins in the N4 size class kg/year
c5<- (6.5*24*365*a)*0.001  #Consumption of urchins in the N5 size class kg/year
c6<- (7.5*24*365*a)*0.001  #Consumption of urchins in the N6 size class kg/year

tset <- seq(from = 0, to = 100, length.out = 100) #time steps



du.simu <- NaN*tset; du.simu[1] <- 0.2 #empty vector for du simulation and initial value for kelp consumption
Kt.simu <- NaN*tset; Kt.simu[1] <- 4000 #empty vector for Kt simulation and initial value for kelp carrying capcity


k.simu <- NaN*tset; k.simu[1] <- 4000 #empty vector for k simulation and initial value for kelp
u1.simu <- NaN*tset; u1.simu[1] <- 20 #empty vector for u1 simulation and initial value of N1 urchins
u2.simu <- NaN*tset; u2.simu[1] <- 15 #empty vector for u2 simulation and initial value of N2 urchins
u3.simu <- NaN*tset; u3.simu[1] <- 30 #empty vector for u3 simulation and initial value of N3 urchins
u4.simu <- NaN*tset; u4.simu[1] <- 15 #empty vector for u4 simulation and initial value of N4 urchins
u5.simu <- NaN*tset; u5.simu[1] <- 15 #empty vector for u5 simulation and initial value of N5 urchins
u6.simu <- NaN*tset; u6.simu[1] <- 5 #empty vector for u6 simulation and initial value of N6 urchins


g1.simu <- NaN*tset; g1.simu[1] <- .99 #empty vector for gamma simulation and initial value for gamma 1
g2.simu <- NaN*tset; g2.simu[1] <- .9 #empty vector for gamma simulation and initial value for gamma 2
g3.simu <- NaN*tset; g3.simu[1] <- .7 #empty vector for gamma simulation and initial value for gamma 3
g4.simu <- NaN*tset; g4.simu[1] <- .5 #empty vector for gamma simulation and initial value for gamma 4
g5.simu <- NaN*tset; g5.simu[1] <- .3 #empty vector for gamma simulation and initial value for gamma 5
g6.simu <- NaN*tset; g6.simu[1] <- .1 #empty vector for gamma simulation and initial value for gamma 6




############### MODEL with No fishing 100X2 M Transect##############################################################################
###################################################################################################################

Z_1<-.9 #annual probability of death
Z_2<-.5 #annual probability of death
Z_3<-.3 #annual probability of death
Z_4<-.3 #annual probability of death
Z_5<-.1 #annual probability of death
Z_6<-.1 #annual probability of death


F_4<- 0 # harvest of N4 urchins every year from kelp site
F_5<-0 # harvest of N5 urchins every year from kelp site
F_6<- 0 # harvest of N6 urchins every year from kelp site



for(i in 2:100){
  
  # calculate change in time
  dt <- tset[i] - tset[i-1]
  
  # store dummy variables

  du<-du.simu[i-1]
  Kt<-Kt.simu[i-1]
  
  k <- k.simu[i-1]
  u1 <- u1.simu[i-1]
  u2 <- u2.simu[i-1]
  u3 <- u3.simu[i-1]
  u4 <- u4.simu[i-1]
  u5 <- u5.simu[i-1]
  u6 <- u6.simu[i-1]
  
  
  gamma_1<- g1.simu[i-1]
  gamma_2<- g2.simu[i-1]
  gamma_3<- g3.simu[i-1]
  gamma_4<- g4.simu[i-1]
  gamma_5<- g5.simu[i-1]
  gamma_6<- g6.simu[i-1]
  

  
  # calculate change in population size for each state variable 
  
  Kt<- 4000 +500*cos(pi*2*tset[i])*dt # seasonal carrying capacity of kelp 
  
  du<-(c1+0.0032674)*u1+
    (c2+0.0032674)*u2+
    (c3+0.0032674)*u3+
    (c4+0.0032674)*u4+
    (c5+0.0032674)*u5+
    (c6+0.0032674)*u6 # amount of kelp urchins are eating/ kelp death 'u'
  
  
  dk<- (((1-(k/Kt)))*k)-du # kelp density
  
  
  gamma_1<- i_gamma_1 
  gamma_2<- i_gamma_2 
  gamma_3<- i_gamma_3 
  gamma_4<- i_gamma_4 
  gamma_5<- i_gamma_5 
  gamma_6<- i_gamma_6 
  
  
  du1 <- (((Ro*((k/K)+.3))*(u4+u5+u6))-(gamma_1*exp(-Z_1)*u1))*dt  
  du2 <-  ((gamma_1*exp(-Z_1)*u1)+(((1-gamma_2)*exp(-Z_2)*u2)))*dt
  du3 <-  ((gamma_2*exp(-Z_2)*u2)+(((1-gamma_3)*exp(-Z_3)*u3)))*dt
  du4 <-  ((gamma_3*exp(-Z_3)*u3)+(((1-gamma_4)*exp(-Z_4)*u4)))*dt
  du5 <-  ((gamma_4*exp(-Z_4)*u4)+(((1-gamma_5)*exp(-Z_5)*u5)))*dt
  du6 <-  ((gamma_5*exp(-Z_5)*u5)+(((1-gamma_6)*exp(-Z_6)*u6)))*dt
  
  
  # calculate total population size
  
  u1.simu[i]<-pmax(0, du1)
  u2.simu[i]<-pmax(0, du2)
  u3.simu[i]<-pmax(0, du3)
  u4.simu[i]<-pmax(0, (du4- F_4))
  u5.simu[i]<-pmax(0, (du5- F_5))
  u6.simu[i]<-pmax(0, (du6- F_6))
  k.simu[i]<-max(0, k.simu[i-1]+ dk)
  du.simu [i]<-du
  Kt.simu [i]<- Kt
  
  g1.simu[i] <- gamma_1
  g2.simu[i]<- gamma_2
  g3.simu[i] <- gamma_3
  g4.simu[i] <- gamma_4
  g5.simu[i] <- gamma_5
  g6.simu[i] <- gamma_6
  
  
  
}

#Plot kelp biomass through time 

plot(x = tset, y = k.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,5000))

#Plot size-specific urchin abundance through time 

plot(x = tset, y = u1.simu, col = '#43225c', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,100),
     ylim=c(0,2000))
lines(x = tset, y = u2.simu, col = '#897295', lwd = 4)
lines(x = tset, y = u3.simu, col = '#8284bc', lwd = 4)
lines(x = tset, y = u4.simu, col = '#cecee9', lwd = 4)
lines(x = tset, y = u5.simu, col = '#695973', lwd = 4)
lines(x = tset, y = u6.simu, col = '#c7b6c8', lwd = 4)
legend(x = max(2)*0.849, y = max(500), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('#43225c', '#897295', '#8284bc','#cecee9', '#695973', '#c7b6c8' ), horiz=TRUE, bty = "n")


library(dplyr)

library(tidyverse)

## MODEL with mild fishing   ####################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################################

k.simu_2 <- NaN*tset; k.simu_2[1] <- 4000



F_4_2<- 1 # harvest of N4 urchins every year from kelp site
F_5_2<-1 # harvest of N5 urchins every year from kelp site
F_6_2<- 1 # harvest of N6 urchins every year from kelp site



for(i in 2:100){
  
  # calculate change in time
  dt <- tset[i] - tset[i-1]
  
  # store dummy variables
  
  du<-du.simu[i-1]
  Kt<-Kt.simu[i-1]
  
  k <- k.simu[i-1]
  u1 <- u1.simu[i-1]
  u2 <- u2.simu[i-1]
  u3 <- u3.simu[i-1]
  u4 <- u4.simu[i-1]
  u5 <- u5.simu[i-1]
  u6 <- u6.simu[i-1]
  
  
  gamma_1<- g1.simu[i-1]
  gamma_2<- g2.simu[i-1]
  gamma_3<- g3.simu[i-1]
  gamma_4<- g4.simu[i-1]
  gamma_5<- g5.simu[i-1]
  gamma_6<- g6.simu[i-1]
  
  
  
  # calculate change in population size for each state variable 
  
  Kt<- 4000 +500*cos(pi*2*tset[i])*dt # seasonal carrying capacity of kelp 
  
  du<-(c1+0.0032674)*u1+
    (c2+0.0032674)*u2+
    (c3+0.0032674)*u3+
    (c4+0.0032674)*u4+
    (c5+0.0032674)*u5+
    (c6+0.0032674)*u6 # amount of kelp urchins are eating/ kelp death 'u'
  
  
  dk<- (((1-(k/Kt)))*k)-du # kelp density
  
  
  gamma_1<- i_gamma_1 
  gamma_2<- i_gamma_2 
  gamma_3<- i_gamma_3 
  gamma_4<- i_gamma_4 
  gamma_5<- i_gamma_5 
  gamma_6<- i_gamma_6 
  
  
  du1 <- (((Ro*((k/K)+.3))*(u4+u5+u6))-(gamma_1*exp(-Z_1)*u1))*dt  
  du2 <-  ((gamma_1*exp(-Z_1)*u1)+(((1-gamma_2)*exp(-Z_2)*u2)))*dt
  du3 <-  ((gamma_2*exp(-Z_2)*u2)+(((1-gamma_3)*exp(-Z_3)*u3)))*dt
  du4 <-  ((gamma_3*exp(-Z_3)*u3)+(((1-gamma_4)*exp(-Z_4)*u4)))*dt
  du5 <-  ((gamma_4*exp(-Z_4)*u4)+(((1-gamma_5)*exp(-Z_5)*u5)))*dt
  du6 <-  ((gamma_5*exp(-Z_5)*u5)+(((1-gamma_6)*exp(-Z_6)*u6)))*dt
  
  
  # calculate total population size
  
  u1.simu[i]<-pmax(0, du1)
  u2.simu[i]<-pmax(0, du2)
  u3.simu[i]<-pmax(0, du3)
  u4.simu[i]<-pmax(0, (du4- F_4_2))
  u5.simu[i]<-pmax(0, (du5- F_5_2))
  u6.simu[i]<-pmax(0, (du6- F_6_2))
  k.simu_2[i]<-max(0, k.simu_2[i-1]+ dk)
  du.simu [i]<-du
  Kt.simu [i]<- Kt
  
  g1.simu[i] <- gamma_1
  g2.simu[i]<- gamma_2
  g3.simu[i] <- gamma_3
  g4.simu[i] <- gamma_4
  g5.simu[i] <- gamma_5
  g6.simu[i] <- gamma_6
  
  
  
}


#Plot kelp biomass through time 

plot(x = tset, y = k.simu_2, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,2000))


#Plot size-specific urhcin abundance through time

plot(x = tset, y = u1.simu, col = '#43225c', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,100),
     ylim=c(0,2000))
lines(x = tset, y = u2.simu, col = '#897295', lwd = 4)
lines(x = tset, y = u3.simu, col = '#8284bc', lwd = 4)
lines(x = tset, y = u4.simu, col = '#cecee9', lwd = 4)
lines(x = tset, y = u5.simu, col = '#695973', lwd = 4)
lines(x = tset, y = u6.simu, col = '#c7b6c8', lwd = 4)
legend(x = max(2)*0.849, y = max(500), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('#43225c', '#897295', '#8284bc','#cecee9', '#695973', '#c7b6c8' ), horiz=TRUE, bty = "n")




## MODEL with heavy fishing   #################################################################################################################################################################################################
k.simu_3 <- NaN*tset; k.simu_3[1] <- 100

F_4_3<- 2
F_5_3<-2
F_6_3<-2



for(i in 2:100){
  
  # calculate change in time
  dt <- tset[i] - tset[i-1]
  
  # store dummy variables
  
  du<-du.simu[i-1]
  Kt<-Kt.simu[i-1]
  
  k <- k.simu[i-1]
  u1 <- u1.simu[i-1]
  u2 <- u2.simu[i-1]
  u3 <- u3.simu[i-1]
  u4 <- u4.simu[i-1]
  u5 <- u5.simu[i-1]
  u6 <- u6.simu[i-1]
  
  
  gamma_1<- g1.simu[i-1]
  gamma_2<- g2.simu[i-1]
  gamma_3<- g3.simu[i-1]
  gamma_4<- g4.simu[i-1]
  gamma_5<- g5.simu[i-1]
  gamma_6<- g6.simu[i-1]
  
  
  
  # calculate change in population size for each state variable 
  
  Kt<- 4000 +500*cos(pi*2*tset[i])*dt # seasonal carrying capacity of kelp 
  
  du<-(c1+0.0032674)*u1+
    (c2+0.0032674)*u2+
    (c3+0.0032674)*u3+
    (c4+0.0032674)*u4+
    (c5+0.0032674)*u5+
    (c6+0.0032674)*u6 # amount of kelp urchins are eating/ kelp death 'u'
  
  
  dk<- (((1-(k/Kt)))*k)-du # kelp density
  
  
  gamma_1<- i_gamma_1 
  gamma_2<- i_gamma_2 
  gamma_3<- i_gamma_3 
  gamma_4<- i_gamma_4 
  gamma_5<- i_gamma_5 
  gamma_6<- i_gamma_6 
  
  
  du1 <- (((Ro*((k/K)+.3))*(u4+u5+u6))-(gamma_1*exp(-Z_1)*u1))*dt  
  du2 <-  ((gamma_1*exp(-Z_1)*u1)+(((1-gamma_2)*exp(-Z_2)*u2)))*dt
  du3 <-  ((gamma_2*exp(-Z_2)*u2)+(((1-gamma_3)*exp(-Z_3)*u3)))*dt
  du4 <-  ((gamma_3*exp(-Z_3)*u3)+(((1-gamma_4)*exp(-Z_4)*u4)))*dt
  du5 <-  ((gamma_4*exp(-Z_4)*u4)+(((1-gamma_5)*exp(-Z_5)*u5)))*dt
  du6 <-  ((gamma_5*exp(-Z_5)*u5)+(((1-gamma_6)*exp(-Z_6)*u6)))*dt
  
  
  # calculate total population size
  
  u1.simu[i]<-pmax(0, du1)
  u2.simu[i]<-pmax(0, du2)
  u3.simu[i]<-pmax(0, du3)
  u4.simu[i]<-pmax(0, (du4- F_4_3))
  u5.simu[i]<-pmax(0, (du5- F_5_3))
  u6.simu[i]<-pmax(0, (du6- F_6_3))
  k.simu_3[i]<-max(0, k.simu_3[i-1]+ dk)
  du.simu [i]<-du
  Kt.simu [i]<- Kt
  
  g1.simu[i] <- gamma_1
  g2.simu[i]<- gamma_2
  g3.simu[i] <- gamma_3
  g4.simu[i] <- gamma_4
  g5.simu[i] <- gamma_5
  g6.simu[i] <- gamma_6
  
}

#Plot kelp biomass through time 

plot(x = tset, y = k.simu_3, col = 'dark green', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,5000))

#Plot size-specific urhcin abundance through time

plot(x = tset, y = u1.simu, col = '#43225c', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,100),
     ylim=c(0,500))
lines(x = tset, y = u2.simu, col = '#897295', lwd = 4)
lines(x = tset, y = u3.simu, col = '#8284bc', lwd = 4)
lines(x = tset, y = u4.simu, col = '#cecee9', lwd = 4)
lines(x = tset, y = u5.simu, col = '#695973', lwd = 4)
lines(x = tset, y = u6.simu, col = '#c7b6c8', lwd = 4)
legend(x = max(2)*0.849, y = max(500), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('#43225c', '#897295', '#8284bc','#cecee9', '#695973', '#c7b6c8' ), horiz=TRUE, bty = "n")






##################################### Plot with kelp biomass at each level of fishing pressure 



plot(x = tset, y = k.simu, col = '#006a4e', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='kelp abundance (kg)',
     xlim=c(0,100),
     ylim=c(0,600))
lines(x = tset, y = k.simu_2, col = '#5ca08e', lwd = 4)
lines(x = tset, y = k.simu_3, col = '#b8d5cd', lwd = 4)
legend(x = max(100)*0.7, y = max(610), legend = c('No fishing','Fishing', 'High Fishing'),lwd = 4, col =c('#006a4e', '#5ca08e', '#b8d5cd'))

view(k.simu_3)
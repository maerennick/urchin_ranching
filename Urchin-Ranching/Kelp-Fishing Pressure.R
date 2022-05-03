##### Different amounts of fishing pressure 



###############Parameters 

K<-1200 #Max kelp 
R<-.5 #growth rate of kelp 
u<-0 #harvest of kelp 
K_u<-600 #max urchins
Ro<-6 #growth rate 
a_1<-1.25
a_2<-0.05
b_1<-.15
b_2<-0.05
alpha<-1


gamma_1<-.99 #100% of 1cm urchins will transition into the 2cm size class within a year
gamma_2<-.9
gamma_3<-0.7
gamma_4<-0.5
gamma_5<-0.3
gamma_6<-0.1


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


k.simu <- NaN*tset; k.simu[1] <- 500
u1.simu <- NaN*tset; u1.simu[1] <- 20
u2.simu <- NaN*tset; u2.simu[1] <- 15
u3.simu <- NaN*tset; u3.simu[1] <- 30
u4.simu <- NaN*tset; u4.simu[1] <- 15
u5.simu <- NaN*tset; u5.simu[1] <- 10
u6.simu <- NaN*tset; u6.simu[1] <- 5


###############No fishing #################################################################################################################################################################################################

Z_1<-.9
Z_2<-.5
Z_3<-.3
Z_4<-.3
Z_5<-.1
Z_6<-.1





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
  k.simu[i]<-k.simu[i-1]+ dk
  du.simu [i]<-du
  dr.simu [i]<-dr
  dd.simu [i]<-dd
  
}

plot(x = tset, y = k.simu, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time (years)',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,2000))


plot(x = tset, y = u1.simu, col = '#43225c', type='l', lwd = 4, las = 1, xlab = 'Time (years)',ylab='Urchin Abundance (no. of urchins)', 
     xlim=c(0,50),
     ylim=c(0,800))
lines(x = tset, y = u2.simu, col = '#897295', lwd = 4)
lines(x = tset, y = u3.simu, col = '#8284bc', lwd = 4)
lines(x = tset, y = u4.simu, col = '#cecee9', lwd = 4)
lines(x = tset, y = u5.simu, col = '#695973', lwd = 4)
lines(x = tset, y = u6.simu, col = '#c7b6c8', lwd = 4)
legend(x = max(2)*0.849, y = max(800), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('#43225c', '#897295', '#8284bc','#cecee9', '#695973', '#c7b6c8' ), horiz=TRUE, bty = "n")






## fishing#################################################################################################################################################################################################

k.simu_2 <- NaN*tset; k.simu_2[1] <- 500


Z_1_2<-.9
Z_2_2<-.5
Z_3_2<-.3
Z_4_2<-(.3+0.4)
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
  
  
  du1 <- ((Ro*(u4+u5+u6)*(1-((u4+u5+u6)/K_u)))-(gamma_1*exp(-Z_1_2)*u1))*dt
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

plot(x = tset, y = k.simu_2, col = 'black', type='l', lwd = 2, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,100),
     ylim=c(0,2000))


plot(x = tset, y = u1.simu, col = '#43225c', type='l', lwd = 4, las = 1, xlab = 'Time (years)',ylab='Urchin Abundance (no. of urchins)', 
     xlim=c(0,50),
     ylim=c(0,800))
lines(x = tset, y = u2.simu, col = '#897295', lwd = 4)
lines(x = tset, y = u3.simu, col = '#8284bc', lwd = 4)
lines(x = tset, y = u4.simu, col = '#cecee9', lwd = 4)
lines(x = tset, y = u5.simu, col = '#695973', lwd = 4)
lines(x = tset, y = u6.simu, col = '#c7b6c8', lwd = 4)
legend(x = max(2)*0.849, y = max(800), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('#43225c', '#897295', '#8284bc','#cecee9', '#695973', '#c7b6c8' ), horiz=TRUE, bty = "n")



#### High Fishing 

k.simu_3 <- NaN*tset; k.simu_3[1] <- 500

Z_1_3<-.9
Z_2_3<-.5
Z_3_3<-.3
Z_4_3<-(.3+0.6)
Z_5_3<-(.1+0.8) #(death rate + fishing pressure)
Z_6_3<-(.1+0.8) #(death rate + fishing pressure)



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
  
  
  du1 <- ((Ro*(u4+u5+u6)*(1-((u4+u5+u6)/K_u)))-(gamma_1*exp(-Z_1_3)*u1))*dt
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

plot(x = tset, y = k.simu_3, col = 'dark green', type='l', lwd = 4, las = 1, xlab = 'Time',ylab='kelp abundance (kg)', 
     xlim=c(0,50),
     ylim=c(0,2000))


plot(x = tset, y = u1.simu, col = '#43225c', type='l', lwd = 4, las = 1, xlab = 'Time (years)',ylab='Urchin Abundance (no. of urchins)', 
     xlim=c(0,100),
     ylim=c(0,800))
lines(x = tset, y = u2.simu, col = '#897295', lwd = 4)
lines(x = tset, y = u3.simu, col = '#8284bc', lwd = 4)
lines(x = tset, y = u4.simu, col = '#cecee9', lwd = 4)
lines(x = tset, y = u5.simu, col = '#695973', lwd = 4)
lines(x = tset, y = u6.simu, col = '#c7b6c8', lwd = 4)
legend(x = max(2)*0.849, y = max(800), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('#43225c', '#897295', '#8284bc','#cecee9', '#695973', '#c7b6c8' ), horiz=TRUE, bty = "n")










plot(x = tset, y = k.simu, col = '#006a4e', type='l', lwd = 4, las = 1, xlab = 'Time (years)',ylab='Kelp Biomass (kg)',
     xlim=c(0,100),
     ylim=c(0,2000))
lines(x = tset, y = k.simu_2, col = '#5ca08e', lwd = 4)
lines(x = tset, y = k.simu_3, col = '#b8d5cd', lwd = 4)
legend(x = max(50)*0.665, y = max(2000), legend = c('No fishing','Mild Fishing', 'Heavy Fishing'),lwd = 4, col =c('#006a4e', '#5ca08e', '#b8d5cd'))


tail(k.simu_3)



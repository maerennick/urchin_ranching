
#### Trying a new appraoch to the model 

gamma_1<-0.99 #100% of 1cm urchins will transition into the 2cm size class within a year
gamma_2<-0.8
gamma_3<-0.7
gamma_4<-0.6
gamma_5<-0.5
gamma_6<-0.4

#Death -> add a death function here as well 

Z_1<-.5
Z_2<-.5
Z_3<-.3
Z_4<-.3
Z_5<-.1
Z_6<-.1


#other parameters

K<-10000
R<-0.5
u<-0
r<-2 #(recruitment function)

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
u2.simu <- NaN*tset; u2.simu[1] <- 100
u3.simu <- NaN*tset; u3.simu[1] <- 100
u4.simu <- NaN*tset; u4.simu[1] <- 100
u5.simu <- NaN*tset; u5.simu[1] <- 100
u6.simu <- NaN*tset; u6.simu[1] <- 100


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
  dk <- ((1-R*(1-(k/K)))*k-u*k-(u1*c1+u2*c2+u3*c3+u4*c4+u5*c5+u6*c6))*dt
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

}


#u1.simu
u2.simu
#u3.simu
#u4.simu
#u5.simu
#u6.simu

plot(x = tset, y = u1.simu, col = 'black', type='p', lwd = 2, las = 1, xlab = 'Time',ylab='urchin abundance', 
     xlim=c(0,30),
     ylim=c(0,500))
lines(x = tset, y = u2.simu, col = 'blue', lwd = 2)
lines(x = tset, y = u3.simu, col = 'green', lwd = 2,lty = 2)
lines(x = tset, y = u4.simu, col = 'dark blue', lwd = 2,lty = 2)
lines(x = tset, y = u5.simu, col = 'light blue', lwd = 2,lty = 2)
lines(x = tset, y = u6.simu, col = 'dark green', lwd = 2,lty = 2)
legend(x = max(10)*0.849, y = max(500), legend = c('u1','u2', 'u3','u4','u5', 'u6'),lwd = 2, col =c('black', 'blue', 'green','dark blue', 'light blue', 'dark green' ))


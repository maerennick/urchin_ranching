#### Urchin Population Model 




### Parameters

r<-0.5 #reproduction rate
Ma1<- 0.01 #mortality rate
Fi<-0.4 #fishing incidental death
Fa<-0.03 #fishing
s<-0.05 # Conversion efficiency of P into D
Din<- 63.38
k<- .327 #growth rate constant
a1<-1
c<-.003482 #size-dependent consumption rate 
g<-.327

  

## delte these later 

P<-2000
D2<- 20
D3<- 20
D4<-20
D5<-20


tset1 <- seq(from=0, to=1000, length.out = 100000)
D.simu1 <- NaN*tset1; D.simu1[1] <- 5

#create holding vectors and set intital conditions 
for(i in 2:length(tset1)){
  dt <- tset1[i]-tset1[i-1]
  
  #store dummy variables
  D1 <- D.simu1[i-1]

  
  #calculate the change in resource and N and M at each timestep
  #dD <- ((r*D2*D3*D4)+ D1*exp(-(Ma1+Fi+Fa))+ s*P*c*D1*(Din*(1-exp(-k*a1)))-g*D1)*dt
  dD<-100
 
  #calcualte the total size of resource N and M
  D.simu1[i] <- D1 + dD
  

}

D.simu1

plot(tset1,D.simu1,type='l',las=1,lwd=2,ylab='Population Size',xlab='Time')



tset2 <- seq(from=0, to=1000, length.out = 100000)
D.simu2 <- NaN*tset1; D.simu2[1] <- 300

#create holding vectors and set intital conditions 
for(i in 2:length(tset2)){
  
  dt <- tset2[i]-tset2[i-1]
  
  #store dummy variables
  D1 <- D.simu2[i-1]
  
  
  #calculate the change in resource and N and M at each timestep
  dD <- (1000*exp(-.22))*dt
  
  #calcualte the total size of resource N and M
  D.simu2[i] <- D1 + dD
  
}

plot(tset2,D.simu2,type='l',las=1,lwd=2,ylab='Population Size',xlab='Time')

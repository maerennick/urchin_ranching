################Random Forrest


rm(list=ls())
require(deSolve)      #ODE solver
require(randomForest) #GSA package
require(rpart)        #Classification/regression tree package
require(rpart.plot)   #makes nicer trees


a_1<-1.25
a_2<-0.05
b_1<-.15
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


threestage <- function (t, x, parms) {
  du = x[1]
  k= x[2]
  u1 = x[3]
  u2 = x[4]
  u3 = x[5]
  u4 = x[6]
  u5 = x[7]
  u6 = x[8]
  
  
  K= parms[1]  
  Ro=parms[2]      
  K_u= parms[3]    
  gamma_1=parms[4]  
  Z_1=parms[5]      
  gamma_2= parms[6] 
  Z_2= parms[7]    
  gamma_3=parms[8]  
  Z_3= parms[9]     
  gamma_4= parms[10]
  Z_4= parms[11]
  Z_5= parms[12]    
  gamma_5= parms[13]
  Z_6= parms[14]    
  gamma_6= parms[15]
  dr= parms[16]     
  dd= parms[17]
  
  
  du<-(c1+0.0032674)*u1+
    (c2+0.0032674)*u2+
    (c3+0.0032674)*u3+
    (c4+0.0032674)*u4+
    (c5+0.0032674)*u5+
    (c6+0.0032674)*u6
  
  dk<- (((1-(k/K))*(dr-dd))*k)-du
  du1 <- ((Ro*(u4+u5+u6)*(1-((u4+u5+u6)/K_u)))-(gamma_1*exp(-Z_1)*u1))
  du2 <-  ((gamma_1*exp(-Z_1)*u1)+(((1-gamma_2)*exp(-Z_2)*u2)))
  du3 <-  ((gamma_2*exp(-Z_2)*u2)+(((1-gamma_3)*exp(-Z_3)*u3)))
  du4 <-  ((gamma_3*exp(-Z_3)*u3)+(((1-gamma_4)*exp(-Z_4)*u4)))
  du5 <-  ((gamma_4*exp(-Z_4)*u4)+(((1-gamma_5)*exp(-Z_5)*u5)))
  du6 <-  ((gamma_5*exp(-Z_5)*u5)+(((1-gamma_6)*exp(-Z_6)*u6)))
  
  return(list(c(du, dk, du1, du2, du3, du4, du5, du6)))
}






sims=4e3
parameters <- data.frame(rbind(
  runif(sims,400,410), #K
  runif(sims,4,5),   #Ro
  runif(sims,400,410),   #K_u
  runif(sims,0.6,0.99),   #gamma_1
  runif(sims,0.8,0.99),   #Z_1
  runif(sims,0.6,0.99), #gamma_2
  runif(sims,0.6,0.99), #Z_2
  runif(sims,0.5,0.99), #gamma_3
  runif(sims,0.6,0.99),    #Z_3
  runif(sims,0.3,0.99),   #gamma_4
  runif(sims,0.6,0.99),   #Z_4
  runif(sims,0.6,0.99),    #Z_5
  runif(sims,0.3,0.99),   #gamma_5
  runif(sims,0.6,0.99),    #Z_6
  runif(sims,0.3,0.99),  #gamma_6
  runif(sims,1.2,1.25),    #dr
  runif(sims,0.1,.15)))  #dd



tf= 200          #run time
times=1:tf

N <-matrix(NA,8,ncol(parameters))  #create vector with NA's to store equilibrium state var. population size
state<-c(du=200, k= 400, u1=20, u2=15, u3=30, u4=15, u5=10, u6=5)     #Kelp forest inits

for(j in 1:ncol(parameters)){
  
state = state*(state>=0) # no negative population sizes allowed
 
  K= parameters[1,j]  
  Ro=parameters[2,j]      
  K_u= parameters[3,j]    
  gamma_1=parameters[4,j]     
  Z_1=parameters[5,j]      
  gamma_2= parameters[6,j]    
  Z_2= parameters[7,j]    
  gamma_3=parameters[8,j]     
  Z_3= parameters[9,j]      
  gamma_4= parameters[10,j]
  Z_4= parameters[11,j]
  Z_5= parameters[12,j]      
  gamma_5= parameters[13,j]
  Z_6= parameters[14,j]      
  gamma_6= parameters[15,j]
  dr= parameters[16,j]      
  dd= parameters[17,j]
  
  
  parms=c(K=K,Ro=Ro, K_u=K_u, gamma_1=gamma_1, Z_1=Z_1,gamma_2=gamma_2, Z_2=Z_2,gamma_3=gamma_3, Z_3=Z_3,gamma_4=gamma_4, Z_4=Z_4,gamma_5=gamma_5, Z_5=Z_5,gamma_6=gamma_6, Z_6=Z_6, dr=dr, dd=dd)  
  
  out=lsoda(state, times, threestage, parms)#, method="lsode",maxsteps=1e4, verbose=TRUE)
  

  N[1,j]=out[tf,2]   #kelp pop matrix
  N[2,j]=out[tf,3]   #U1 pop matrix
  N[3,j]=out[tf,4]   #U2 pop matrix
  N[4,j]=out[tf,5]   #U3 pop matrix
  N[5,j]=out[tf,6]   #U4 pop matrix
  N[6,j]=out[tf,7]   #U5 pop matrix
  N[7,j]=out[tf,8]   #U6 pop matrix
  
}


##Set up data for kelp GSA
kelpGSA<-matrix(NA,sims,29)   #dummy matrix to hold state variables and sampled parameter values
kelpGSA[,1:5]<-t(N)          #input state variable values from simulations
kelpGSA[,6:29]<-t(parameters) #input randomly chosen parameter values
kelpGSA<-as.data.frame(kelpGSA)
colnames(kelpGSA)<-c("kelp", "urchin 1","urchin 2","urchin 3","urchin 4","urchin 5","urchin 6", "alpha 1","alpha 2","beta 1","beta 2","Kelpcc","dr","dd","du","Ro","Urchincc","gamma1", "death1","gamma2", "death2","gamma3", "death3","gamma4", "death4","gamma5", "death5","gamma6", "death6")

##Calculate MTL for kelp
MTLKelp<- matrix(NA,1,nrow(kelpGSA))  #dummy df for mean trophic level output
MTLcalcdata<-as.data.frame(t(kelpGSA[,1:5]))
MTLcalcdata<-MTLcalcdata[-1,]  #drop kelp biomass from pop size df
MTLcalcdata<-rbind(MTLcalcdata,colSums(MTLcalcdata)) #add total biomass of inverts as new row
## now calculate mean trophic level
MTLKelp<- (2*MTLcalcdata[1,]/MTLcalcdata[5,]) +(2*MTLcalcdata[2,]/MTLcalcdata[5,])+
  (2*MTLcalcdata[3,]/MTLcalcdata[5,]) + (2.5*MTLcalcdata[4,]/MTLcalcdata[5,])
MTLKelp<-t(MTLKelp)  #need to transpose MTL to match the structure of the GSA parameter object

#Run Kelp GSA
MTLkelpGSA<-cbind(MTLKelp,kelpGSA) #combine MTL values with parameter values
colnames(MTLkelpGSA)[1] <- "MTL"

##Run random forest on mean trophic level data
mtlkelp.rf<-randomForest(MTLkelpGSA[,1] ~ .,
                         data       = MTLkelpGSA[,7:30],
                         importance = TRUE)
print(mtlkelp.rf)
round(importance(mtlkelp.rf),2)
mtlkelp.imp = pmax (importance (mtlkelp.rf,scale = FALSE) [, 1], 0)
mtlkelp.imp = mtlkelp.imp/sum(mtlkelp.imp)


#Bar plot of normalized importance values
pdf("ImportanceMTLKelp.pdf")
par(mai=c(1.02,1.02,0.62,0.42))
barplot (mtlkelp.imp,
         horiz = T,
         main  = 'Mean trophic level',
         xlim  = c (0, .5),  las=1)
dev.off()























###################################
### Bifurcation Diagram 
Kset <- seq(from = 0.01, to = 1, length.out = 200)


u1starset<- rep(NaN, times = length(Kset))
u2starset<- rep(NaN, times = length(Kset))
u3starset<- rep(NaN, times = length(Kset))
u4starset<- rep(NaN, times = length(Kset))
u5starset<- rep(NaN, times = length(Kset))
u6starset<- rep(NaN, times = length(Kset))
kstarset<- rep(NaN, times = length(Kset))
dustarset<- rep(NaN, times = length(Kset))
drstarset<- rep(NaN, times = length(Kset))
ddstarset<- rep(NaN, times = length(Kset))
utotalstarset<- rep(NaN, times = length(Kset))


timepoints <- seq(from = 1, to = 100, length.out = 200)


for(j in 1:length(Kset)) {
  Z_5 <- Kset[j]
  
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
  utotal.simu <- NaN*tset; u6.simu[1] <- 95
  
  
  for(i in 2:length(timepoints)){
    
    dt <- timepoints[i]-timepoints[i-1]
    
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
    utotal <- utotal.simu[i-1]
    
    
    
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
    
    dutotal<-u1+u2+u3+u4+u5+u6
    
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
    utotal.simu [i] <- dutotal
    
  }
  
  
  u1starset[j]<- u1.simu[length(timepoints)]
  u2starset[j]<- u2.simu[length(timepoints)]
  u3starset[j]<- u3.simu[length(timepoints)]
  u4starset[j]<- u4.simu[length(timepoints)]
  u5starset[j]<- u5.simu[length(timepoints)]
  u6starset[j]<- u6.simu[length(timepoints)]
  kstarset[j]<-  k.simu [length(timepoints)]
  dustarset[j]<- du.simu[length(timepoints)]
  drstarset[j]<- dr.simu[length(timepoints)]
  ddstarset[j]<- dd.simu[length(timepoints)]
  utotalstarset[j]<-utotal.simu[length(timepoints)]
  
}


#State R 1 Bifurcation Diagram 

par(mfrow=c(1,3))

plot(x = Kset, y = utotalstarset, main="1", 
     type='l', 
     lwd=2, 
     xlab = "Ro", 
     ylab = "Carrying capacity of u1")
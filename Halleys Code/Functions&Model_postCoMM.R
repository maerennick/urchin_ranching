##***PARAMETERS***  
# Mvec =  vector of motalities from age i to the next (annual mortality)
# eggPf = proportion of female eggs = 0.5
# a = BH density independent mortality
# b = BH density dependent mortality
# Ro = recruitment (age 1) in unfished conditions
# h = steepness
# EPR = Eggs-per-recruit in unfished conditions
# alpha= half-saturation constant for 'effective reproductive ratio' (try values betweem 0.1-0.2)  
# effRatio ='effective reproductive ratio' 
# Fec = fecundity vector at age i = 1.5-2.5million eggs per year
# Mature = proportion of mature age i vector
# delta = probability of incidental death
# Omega = probability of capture retention
# Omax = maximum probability of functional capture retention
# u = vulnerability to capture (= E * q)
# q = TRUE catchability
# q_max = MAX TRUE catchability with hypoxia (when hypoxic depth = 0m)
# E = fishing effort (number of pots per week; com = commercial tribal, rec = state recreational)
# Fr = fishing mortality from retention
# Fi = fishing moatlity from incidental catch
# d = constant for catchability function (q) given hypoxia
#####################################################################


##Quick functions to calculate Nt+1 and Catch, based on input vectors of abundance, mortality and fishing (Barnov instantaneous catch equation)
step.forward<-function(N,M,u,Omega,delta) {
  Fr=u*Omega
  Fi=u*(1-Omega)*(delta)
  return(N*exp(-(M+Fr+Fi)/52))
}

step.handle.mort<-function(N,M,u,Omega,delta) {
  Fr=u*Omega
  Fi=u*(1-Omega)*(delta)
  return(N*Fi/(M+Fr+Fi)*(rep(1,length(N))-exp(-(M+Fr+Fi)*1/52)))
}


step.catch<-function(N,M,u,Omega,delta){
  Fr=u*Omega
  Fi=u*(1-Omega)*(delta)
  return(N*Fr/(M+Fr+Fi)*(rep(1,length(N))-exp(-(M+Fr+Fi)*1/52)))
}


##Catchability function to determine q , with q = 0.005 when hypoxia depth = 0m (i.e., q_max)
q.func<-function(DO,hypoxia,q_max,qmvec){
  ifelse(hypoxia==1 && DO<= 30,return(d*DO + q_max),return(q))
}


##Control illegal retention relationship 
ill.func<-function(ill,Omega,Omax,catch,total.effort,c){
  if(ill==1 && total.effort>0)
    return(Omax*exp(-(catch/total.effort)/c)) 
   else
    return(Omega)
}


##Function for Effective Reproductive Ratio. Bio==1 turns on Eff. Ratio Function
eff.reprd.func<-function(Nf,Fec,Mature,effRatio,bio){
  if(bio==1)
    eggs <- sum(Nf*Fec*effRatio*Mature)
  else
    eggs <- sum(Nf*Fec*Mature) 
  return(eggs)
}



##Dungeness crab model*****
dungy.age.model<-function(Nm,Nf,Mmvec,Mfvec,qmvec,qfvec,Omegamvec,deltamvec,Omegafvec,deltafvec,E.rec,E.com,Max.Catch,Omega,Omax,hypoxia,ill,c) {
  
  # Cycling through fishing weeks
  cum.catch<-0 
  cum.catch.fem<-0
  cum.catch.sub<-0
  cum.handle.fem<-0
  cum.handle.sub<-0
  
  
  for (i in 1:n.weeks){
 
    #Catchability. Hypoxia==1 --> non-constant q
    do.depth<-DO$Depth[i]
    qmvec[3:10]<-q.func(do.depth,hypoxia,q_max,qmvec) 
    qfvec[4:10]<-q.func(do.depth,hypoxia,q_max,qfvec)                                 
    
    #This mimics the tally of of Summer catch that dictates opening the winter season
    if(i>10 && cum.catch>=Max.Catch){
      E.com.2.use = 0
      E.rec.2.use = 0
    }
    else{
      E.com.2.use <- E.com[i]
      E.rec.2.use <- E.rec[i]
    }
    
    E.total=E.com.2.use+E.rec.2.use #for illegal take function
    
    # Here calculate u (probability of capture), based on E (commercial & recreational) and q (true catchability)
    umvec<-(E.com.2.use+E.rec.2.use)*qmvec
    ufvec<-(E.com.2.use+E.rec.2.use)*qfvec
    
    #Males
    Nmt.plus.1<-step.forward(Nm,Mmvec,umvec,Omegamvec,deltamvec)
    
    handle.mort.sub<-sum(step.handle.mort(Nm,Mmvec,umvec,Omegamvec,deltamvec))
    cum.handle.sub<-cum.handle.sub+handle.mort.sub
    
    Nmt.catch<-step.catch(Nm,Mmvec,umvec,Omegamvec,deltamvec)
    cum.catch<-cum.catch+sum(Nmt.catch[4:10])
    cum.catch.sub<-cum.catch.sub+Nmt.catch[3]
    
    #Females
    Nft.plus.1<-step.forward(Nf,Mfvec,ufvec,Omegafvec,deltafvec)
    
    handle.mort.fem<-sum(step.handle.mort(Nf,Mfvec,ufvec,Omegafvec,deltafvec))
    cum.handle.fem<-cum.handle.fem+handle.mort.fem
    
    Nft.catch<-sum(step.catch(Nf,Mfvec,ufvec,Omegafvec,deltafvec))
    cum.catch.fem<-cum.catch.fem+Nft.catch
    
    #total number of handling mortalities
    cum.handle.mort<-cum.handle.sub + cum.handle.fem
    
    Nm<-Nmt.plus.1
    Nf<-Nft.plus.1
    
    #Functional illegal retention
    Omegamvec[3] <- ill.func(ill,Omega,Omax,cum.catch,E.total,c) #males
    #Omegafvec[4:10]<- ill.func(ill,Omega,Omax,cum.catch,E.total,c) #females
  }
  
  # Finish projecting forward to the end of the year
  n.weeks.remaining<-52-n.weeks
  Nm<-Nm*exp(-Mmvec*n.weeks.remaining/52)
  Nf<-Nf*exp(-Mfvec*n.weeks.remaining/52)    
  
  return(list(Nm=Nm,Nf=Nf,Catch=cum.catch,Illegal.fem=cum.catch.fem, Illegal.sub=cum.catch.sub,Handling.mort=cum.handle.mort))
  
}


#dungy.age.model(Nm,Nf,Mmvec,Mfvec,qmvec,qfvec,Omegamvec,deltamvec,Omegafvec,deltafvec,E.rec,E.com,Max.Catch,Omega,Omax,hypoxia,ill,c)


##***Recruitment Function*******************************************************
recruitment.fun<-function(Nm,Nf,alpha=0.1,a,b,Fec,Mature,eggPf=0.5,bio){
  
  eggPm=1-eggPf
  ##Proportion of larger (older) mature males for age i mature females in each age class
  ## E.g., age 2 females can capulate with age2-10, age 3 females with age3-10, etc.
  ## Make maximum value 1, otherwise equal to proportion
  ## Calculate effective reproductive ratio (based on Mondo equation)
  Nm.mx<-as.matrix(Nm)
  Nf.mx<-as.matrix(Nf)
  propMF<-matrix(NA, nrow=10)
  
  for(i in 1:10){
    propMF[i,] <-  sum(Nm.mx[-c(1:i),])/Nf.mx[i,]    
    if(propMF[i,] > 1) propMF[i,] = 1
   } 
  propMF=as.vector(propMF)
  effRatio <- propMF/(rep(alpha,10)+propMF)
  
  eggs=eff.reprd.func(Nf,Fec,Mature,effRatio,bio)   #bio==1 --> effective ratio 'on'
  
  ## Function to calculate recruitment (age 0 to 1) with density dependance using Beverton-Holt model:
  ## num of eggs surviving year to year. Exponentiate v for lognormal error 
  R.f = eggPf*(eggs/((a+b*eggs)))
  R.m = eggPm*(eggs/((a+b*eggs)))
  
  return(list(R.m=R.m,R.f=R.f,eggs=eggs))
}

#recruitment.fun(Nm=Nm,Nf=Nf,alpha=0.1,a=a,b=b,Fec=Fec,eggPf=0.5, Mature=Mature, bio=bio)



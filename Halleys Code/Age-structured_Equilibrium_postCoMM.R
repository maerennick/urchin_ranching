rm(list=ls())  
workdir <-"~/Dropbox/HW/Ch 4/ORCA data"
setwd(workdir) 


##***PARAMETERS***  
# Mvec = probability vector of motality from age i to the next (annual mortality)
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
# c = constant for illegal retention
# d = constant for catchability function (q) given hypoxia
#####################################################################

#Source in the Model and associated functions
source('~/Dropbox/HW/Ch 4/R code/Functions/Functions&Model_postCoMM.R', encoding='UTF-8')

##Get a sense of crabs being removed****
#Total proportion removed of each age class
#E.rec= 3000   # 1800 - 3800pots
#E.com = 4000
#q=c(0,0,rep(0.0074,8))
#umvec<-(E.com+E.rec)*q
#Fr=umvec*Omegamvec
#Fi=umvec*(1-Omegamvec)*(deltamvec)

#Instananeous fishing mortality
#(Fr+Fi)*(1/52)
#Total (natural + fishing) instananeous mortality
#(Mmvec+Fr+Fi)*(1/52)
#Propotion removed
#1-exp(-(Mmvec+Fr+Fi)/52) #All
#1-exp(-(Fr+Fi)/52) #Just fishing removals
#1-exp(-Fi/52) #Incidental removals

#**************Set up parameters of simulation**************************************************


# Starting vectors of Nf, Nm*******     
Nm<-c(3500000,900000,200000,100000,50000,25000,11000,5000,2000,1000)
Nf<-c(3500000,900000,200000,100000,50000,25000,11000,5000,2000,1000)


# Natural mortality
Mmvec<-c(0.8,0.8,rep(0.3,8))
Mfvec<-c(0.8,0.8,rep(0.3,8))

## Catchability  (influences probability of capture (u))
all.DO<-read.csv('hypoxic_depth.csv', header = TRUE)
head(all.DO)
DO<-subset(all.DO, all.DO$Type_num==5)

#constant for q function (see excel "CATCHABILITY")
d = -0.00015
#q = d*depth + q_max --> how q was calculated from min hypoxic depth
q_max = 0.005
#q in the absence of hypoxia
q=0.00035
qmvec<-c(0,0,rep(0.00035,8)) 
qfvec<-c(0,0,0,rep(0.00035,7))


## Probability of retention
  #NO ILLEGAL conditions
   Omega=0
   Omegamvec<-c(0,0,0,1,1,1,1,1,1,1)
   Omegafvec<-rep(0,10)

## Probability of incidental death (handling mort.)
  # NO INCIDENTAL MORTALITY
   deltamvec<-rep(0,10) #No handling mortality
   deltafvec<-rep(0,10) 



# Set up parameters of fishing season (effort = # pots per week)
n.weeks=20
E.rec.max= 3000   #low = 3000 ; high = 16000
E.rec<-rep(E.rec.max,n.weeks)
E.com.max = 4000
E.com<-c(rep(E.com.max,n.weeks/4),rep(0,n.weeks/4),rep(E.com.max,n.weeks/4),rep(0,n.weeks/4))
E.com

Erec.seq = seq(0,16000,150) #for Catch vs. Effort graph (allowing no. rec pots to increase, while tribal comm. maxes out at 4000)
Ecom.seq = seq(0,4000,150)
length.com=length(Erec.seq)-length(Ecom.seq)
constant.Ecom<-as.matrix(rep(4000,length.com))
Ecom.seq=as.matrix(Ecom.seq)
Ecom.seq = rbind(Ecom.seq,constant.Ecom)
Ecom.seq = as.vector(Ecom.seq)


# Set Quotas
Max.Catch<-1000000000000000000 #For equlibrium



##B-H parameters
Mature = c(0,0.2,1,1,1,1,1,1,1,1)
Fec=c(0,1000000,2000000,2000000,2000000,1000000,1000000,1000000,1000000,1000000)

Ro = 7000000  #Number of recruits (age 1) without fishing

  #Calculate EPR
  Prop_of_Nf <- matrix(NA,nrow=1,ncol=10)
  Mfmx<- as.matrix(Mfvec)
  Prop_of_Nf[1] = 0.5
  for(r in 2:10){
    Prop_of_Nf[r] <- Prop_of_Nf[r-1]*exp(-Mfmx[r-1])  
  }
  Prop_of_Nfvec<-as.vector(Prop_of_Nf)
  EPR=round(sum(Prop_of_Nfvec*Mature*Fec))

EPR
h=0.5
a = ((1-h)/(4*h))*EPR
b = (5*h-1)/(4*h*Ro)


#others Parameters
alpha=0.1
eggPf=0.5
Omax=0.09 #low = 0.02
c=10

#************************TURN ON or OFF function in the model**********************
bio=0  # 1 = effective reproductive ratio; 0 = no ERR
hypoxia=1 # 1 = non-linear q; 0 = constant q
ill=0 # 1 = functional illegal take (Omgea); 0 = constant illegal take
#**********************************************************************************


n.years<-50
##Create 'NA' bank to store yearly outputs of Nf, Nm, and catch
N.f.output<-matrix(NA,nrow=n.years+1,ncol=10)
N.m.output<-matrix(NA,nrow=n.years+1,ncol=10)
catch.output<-matrix(NA,nrow=n.years,ncol=1)
illegal.fem.output<-matrix(NA,nrow=n.years,ncol=1)
illegal.sub.output<-matrix(NA,nrow=n.years,ncol=1)
egg.output<-matrix(NA,nrow=n.years,ncol=1)
catch.vs.effort<-matrix(NA,nrow=length(Erec.seq),ncol=2)

#Initial conditions (year 1 form stable-age equlibrium excel spreadsheet)
N.f.output[1,]<-Nf
N.m.output[1,]<-Nm



##Run Model*****************************************************************************
quartz()
#par(mfrow=c(3,2))
par(mar=c(5,7,4,3)+0.1)  #c(bottom, left, top, right) 
#plot(c(),c(), xlim=c(0,20000), ylim=c(0,800000), type='l', xlab='', ylab= '', lwd=2, las=1, main="Vary Omega: Effective Ratio", cex.axis=1.2)
#mtext("Catch (no. crabs)", side=2, line=4)
#mtext("Effort (no. pots per week)", side=1, line=2.5)
#legend("topleft", legend=c("0","0.2","0.4","0.6","0.8"), col='black', lty=1, lwd=2, bty='n')
#legend("topleft", legend=c("h=0.8","h=0.6","h=0.5"), col='black', lty=1, lwd=2, bty='n')
#legend("topright", legend=h, bty='n')

#looping over steepness
#for(h in c(0.8,0.6,0.5)){
#h=h
#a = ((1-h)/(4*h))*EPR
#b = (5*h-1)/(4*h*Ro)


plot(c(),c(), xlim=c(0,20000), ylim=c(0,700000), type='l', xlab='', ylab= '', 
     lwd=3, las=1, main="", cex.axis=1.4)
mtext("Crabbing Effort (no. pots week-1)", side=2, line=5, cex=1.5)
mtext("Legal Catch (no. crabs)", side=1, line=2.5, cex=1.5)
legend("topright", legend="δf", bty='n', cex=1.5) #Base-case, Φ, Ω, δ
legend("topleft", legend="f)", bty='n', cex=2)
legend("bottomright", legend=c("no hypoxia","hypoxia"), col=c("gray","black"), lwd=5,lty=1,bty='n', cex=1.2)

i = 0.8 #use i to control Omega for constant ill., Omax for ill function, and the two delta vec

#Turn on to loop over differ constant values for various parameters (Omega - illegal retention, delta - incidental mort, etc.) 0.8 for parameters
#for(i in seq(from=0, to=0.8, by=0.8)){
  #Omega=i
  #deltamvec<-c(0,0,i,i,i,i,i,i,i,i) 
  deltafvec<-c(0,0,i,i,i,i,i,i,i,i)
  #Omax=i

  
#Turn this on to create catch vs. effort plots (also un-# matrices below to save data for plot)
 for (e in 1:length(Erec.seq)){
   E.rec = rep(Erec.seq[e],n.weeks)
   E.com = c(rep(Ecom.seq[e],n.weeks/4),rep(0,n.weeks/4),rep(Ecom.seq[e],n.weeks/4),rep(0,n.weeks/4))

   #Reset Initial conditions
   Nm<-c(3500000,900000,200000,100000,50000,25000,11000,5000,2000,1000)
   Nf<-c(3500000,900000,200000,100000,50000,25000,11000,5000,2000,1000)

   for (t in 1:n.years){
     
    #First thing, go crabbing
    N.tplus.1<-dungy.age.model(Nm=Nm,Nf=Nf,Mmvec=Mmvec,Mfvec=Mfvec,qmvec=qmvec,qfvec=qfvec,Omegamvec=Omegamvec,Omegafvec=Omegafvec,Omega=Omega,
                             deltamvec=deltamvec,deltafvec=deltafvec,E.rec=E.rec,E.com=E.com,Max.Catch=Max.Catch,Omax=Omax,hypoxia=hypoxia,ill=ill,c)
    catch.output[t]<-N.tplus.1$Catch
    illegal.fem.output[t]<-N.tplus.1$Illegal.fem
    illegal.sub.output[t]<-N.tplus.1$Illegal.sub
    
    #Assign zero to age 1 year t+1 (simply a placeholder for post-fishing age 1 recruits)
    N.f.output[t+1,1]<-0
    N.m.output[t+1,1]<-0
    # Update females in output matrix
    N.f.output[t+1,2:10]<-N.tplus.1$Nf[1:9]
    # Update males in output matrix
    N.m.output[t+1,2:10]<-N.tplus.1$Nm[1:9]
    #Update the Nf and Nm vectors to calculate post-crabbing recruitment
    Nf<-N.f.output[t+1,]
    Nm<-N.m.output[t+1,]
    
    
    #Generate recruitement AFTER fishing
    R.tplus.1<-recruitment.fun(Nm=Nm,Nf=Nf,alpha=alpha,a=a,b=b,Fec=Fec,eggPf=eggPf,Mature=Mature,bio=bio)
    # Now add in recruits, that now become age 1
    N.f.output[t+1,1]<-R.tplus.1$R.f
    N.m.output[t+1,1]<-R.tplus.1$R.m
    #Update the Nf and Nm vectors for the next loop iteration
    Nf<-N.f.output[t+1,]
    Nm<-N.m.output[t+1,]
    #Egg output for reference
    #egg.output[t,1]<-R.tplus.1$eggs
   }

  catch.vs.effort[e,1]<-Erec.seq[e]+Ecom.seq[e]
# catch.vs.effort[e,1]<-Erec.seq[e] + E.com.max # constant commercial throughout
  catch.vs.effort[e,2]<-catch.output[n.years]

 }

lines(catch.vs.effort,col="black", lwd=7, lty=1) #black = hypoxia, gray = no hypoxia
#}

#}

catch.vs.effort


















##*****PLOT TIME!************

#Plot catch vs. effort
max(catch.vs.effort[,2])
windows(7,7)
plot(c(),c(), xlim=c(0,20000), ylim=c(0,600000), type='l', xlab='Effort (no. pots per week)', ylab= 'Catch (no. crabs)', lwd=2)

lines(catch.vs.effort,col="black", lwd=2, lty=1)
lines(catch.vs.effort,col="black", lwd=2, lty=2)
lines(catch.vs.effort,col="blue", lwd=2, lty=1)
lines(catch.vs.effort,col="red", lwd=2, lty=1)
#lines(catch.vs.effort,col="orange", lwd=2, lty=1)

legend("topleft", legend=c("200k quota Handling Mort. 0.01-0.8"), col=c("black","black"), lty=c(1,1,1,1), lwd=2, bty='n')
legend("bottomright", legend=c("Constant Illegal Retention", "Sublegal 1% Handling Mort"), col=c("black","blue"), lty=c(1,1), lwd=2, bty='n')

legend("topright", legend=c("h = 0.8, No Effective Ratio"), bty='n')


#Quad plot to view Population abundance, Pop. growth rate, Catch, and Illegal catch       
Nf.total<-rowSums(N.f.output)
Nm.total<-rowSums(N.m.output)

N.total <- Nf.total+Nm.total
N.total

growth.rate <- rep(NA, n.years-1) 
for(i in 1:n.years-1){
growth.rate[i] <- N.total[i+1]/N.total[i]
}

windows(5,5)
par(mfrow=c(2,3))
plot(N.total, xlab='Years',type='l', ylim=c(0,max(N.total)), ylab="Population Abundance")
plot(growth.rate, xlab='Years', type='l', ylim=c(0,2), ylab="Population Growth Rate")
plot(catch.output,xlab='Years', type='l',ylim=c(0,max(catch.output)), ylab="Catch (no. crabs)", las=1)
plot(illegal.fem.output,xlab='Years', type='l',ylim=c(0,max(illegal.fem.output)), ylab="Illegal Female Catch")
plot(illegal.sub.output,xlab='Years', type='l',ylim=c(0,max(illegal.sub.output)), ylab="Illegal Sublegal Catch")




##Just checking stuff....
# Reference: 1% of potrol records were illegal females
illegal.output/catch.output

#eff
growth.rate

#Equilbrium egg production should be around 3E+11 (w/o fishing)
#egg.output
round(N.m.output)
round(N.f.output)

#Mean max catch with effective ratio
mean(c(521974.8,527031.3,532082.5))
sd(c(521974.8,527031.3,532082.5))


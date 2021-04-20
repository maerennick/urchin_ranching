####### Mae Rennick
###### Urchin Growth Rates

##Kenner 1992 (Population dynamics of the sea urchin Strongylocentrotus purpuratus in a Central California kelp forest: recruitment, mortality, growth, and diet) Growth Curve 

## Von Bertalanffy Equation 
#Test diameter as a product of age

Din<- 63.38
K<- .327 #growth rate constant
a0<-0

diameter <- function(a) {
  D <- Din*(1-exp(-K*(a-a0)))
  return(D)
}

age<-c()

diameter(age)


plot(diameter)

#Dd/Da=K(Din-D)
#seasonally adjusted Bertalanffy exists

##for loop where each age is i a0 can tay at zero 
#but then you need to add change in diameter to original diameters unit of time will be one year unless you want to divide it by 12 to have it by month


#size dependent consumption rate based of of mesocosum expiriments

b<-.00089
a<-.003482

#per_cap_consumption<-a+b*D




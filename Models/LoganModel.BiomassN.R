#Creation of a model to determine what significantly impacts the biomassN of the observations. 
require(ggplot2)
require(GGally)
require(viridis)
require(boot)
require(Sleuth3)
require(ResourceSelection)
library(lmtest)
library(glmmTMB)
library(DHARMa)
require(AER)
library(sjstats)
library(MASS)
library(merTools)
library(performance)
library(lmerTest)
library(car)

#data exploration 

#We need to create a dataset with everything but create a new column that makes soil depth into 
#a factor (low, medium, high). 

#Check correlation and remove correlated variables 

#Still may need to standardize scales. 
#One interaction: soil depth with some other variable. 

modified1 <- read.csv("MTNYCData_modified1.csv")

#screen for correlation 
source("screen_cor.R")
colnames(modified1)
screen.cor(modified1[,-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 18, 19, 20, 24, 26, 29, 30, 31, 32, 35)])


#Stepwise reduction for BiomassN

m1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + 
              NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + 
              (1|Site) + (1|Plot) + (1|Year), data = modified1) 
summary(m1)
plot(m1)
vif(m1)
isSingular(m1)

m2 <- lmer(BiomassN ~ factor(diversity) + factor(Season) + Respiration + NO2_NO3 + 
              NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + 
              (1|Site) + (1|Plot) + (1|Year), data = modified1)
summary(m2)

m3 <- update(m2,.~.-RootMass_g)
summary(m3)

m4 <- update(m3,.~.-factor(diversity))
summary(m4)
isSingular(m4)

m5 <- update()

#TESTING 
TEST <- lmer(BiomassN ~ factor(Season) + Respiration + NO2_NO3 + 
               NH4 + Nitrification + Moisture_g + factor(coresection) + 
             (1|Site) + (1|Plot) + (1|Year), data = modified1)

#finding what variables are causing singularity. 

# calculate the matrix rank and determinant
rankm4 <- qr(m4)$rank
detm4 <- det(m4)

# print the rank and determinant
print(paste0("Matrix rank: ", rankX))
print(paste0("Matrix determinant: ", detX))



#m4 is all significant 
plot(m4)

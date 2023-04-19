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

#data exploration 

#We need to create a dataset with everything but create a new column that makes soil depth into 
#a factor (low, medium, high). 

#Check correlation and remove correlated variables 

#Still may need to standardize scales. 
#One interaction: soil depth with some other variable. 

tree.0.10 <- read.csv("tree_0-10.csv")

modified <- read.csv("MTNYCData_modified.csv")

#screen for correlation 
source("screen_cor.R")
colnames(modified)
screen.cor(modified[,-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 23, 25, 28:34)])

#Stepwise reduction for BiomassN

m1 <- lmer(BiomassN ~ success + diversity + Month + Season + Respiration + NO2_NO3 + 
              NH4 + TIN + Mineralization + Nitrification + DEA + RootMass_g + Moisture_g +
              (1|Site) + (1|Plot) + (1|Year), data = tree.0.10) 
summary(m1)

#Remove month from fixed effects and set it as a random effect 

m6 <- lmer(BiomassN ~ NO2_NO3 + NH4 + TIN + Mineralization + Nitrification +  
  DEA + (1 | Site) + (1 | Plot) + (1 | Year) + (1|Month), data = tree.0.10)
summary(m6)


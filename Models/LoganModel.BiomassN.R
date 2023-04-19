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

git config pull.rebase false


#
tree.0.10 <- read.csv("tree_0-10.csv")

#Stepwise reduction for BiomassN

m1 <- lmer(BiomassN ~ success + diversity + Month + Season + Respiration + NO2_NO3 + 
              NH4 + TIN + Mineralization + Nitrification + DEA +
              (1|Site) + (1|Plot) + (1|Year), data = tree.0.10) 
summary(m1)

m2 <- update(m1,.~.-Season)
summary(m2)

m3 <- update(m2,.~.-diversity)
summary(m3)

m4 <- update(m3,.~.-success)
summary(m4)

m5 <- update(m4,.~.-Respiration)
summary(m5)

#Remove month from fixed effects and set it as a random effect 

m6 <- lmer(BiomassN ~ NO2_NO3 + NH4 + TIN + Mineralization + Nitrification +  
  DEA + (1 | Site) + (1 | Plot) + (1 | Year) + (1|Month), data = tree.0.10)
summary(m6)


# this is a testing file

# see if my repository works

library(here)
library(tidyverse)

#data preparation
df <- read.csv("Data/MTNYCData.csv")
df$success <- ifelse(df$Success == "Low", 0, 1)
df$diversity <- ifelse(df$Diversity == "Low", 0, 1)

i <- 1

for (i in 1:392) {
  df$Obs[i] = i
}

m1 <- glm(success ~ Diversity + (1|Obs), data = df, family = binomial)
summary(m1)

library(lme4)
m2 <- lme4::glmer(success~(1|Obs)+
                    BiomassN+BiomassC+RootMass_g+Moisture_g+
                    Nitrification+Mineralization+NO2_NO3, family=binomial, data=df)
summary(m2)

m3 <- lme4::glmer(success~(1|Obs)+as.factor(Season)+
                    BiomassN+BiomassC+RootMass_g+Moisture_g+
                    Nitrification+Mineralization+NO2_NO3, family=binomial, data=df)
summary(m3)


m4 <- lme4::glmer(success~(1|)+
                    BiomassN+BiomassC+RootMass_g+Moisture_g+
                    Nitrification+Mineralization+NO2_NO3, family=binomial, data=df)
summary(m4)

m5 <- lme4::glmer(success~(1|Obs)+as.factor(Season)+
                    BiomassN+BiomassC+RootMass_g+Moisture_g+
                    Nitrification+Mineralization+NO2_NO3, family=binomial, data=df)
summary(m5)

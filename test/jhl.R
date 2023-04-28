# this is a testing file

# see if my repository works

library(here)
library(tidyverse)

#data preparation
df <- read.csv("Data/MTNYCData_modified.csv")
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

## thoughts
# Bidirectional causal effect: C/N -> success (short term); success -> C/N (long term)
# in the short time period restricted by the dataset, it is better to detect C/N -> success channel
# since the C/N concentrations are accumulated in a long period, is less likely to be influenced by
# the one-time plantation activity
# in the discussion, we can talk about the future effect 

df_list <- split(df, df$CoreSection_cm)

library(lme4)
library(lmerTest)
reg_results <- lapply(df_list, function(bioc) {
  lmer(TotalSoilC ~ (1|Site) + success, data = bioc)
})

df.10_30 <- read.csv("Data/tree_10-30.csv")
lm1 <- lmer(TotalSoilC ~ (1|Site) + factor(Month) + factor(Success) + factor(Diversity) + RootMass_g+Moisture_g+Respiration, data = df.)
summary(lm1)

lm2 <- lmer(BiomassC ~ (1|Site) + factor(Month) + factor(Success) + factor(Diversity) + RootMass_g+Moisture_g+Respiration, data = df.10_30)
summary(lm2)

lm(BiomassC ~  factor(success) + factor(Diversity)+ RootMass_g + Moisture_g + Respiration, data = df)

chisq.test(table(df.10_30$TotalSoilC_percent, df.10_30$success))

library(ggplot2)

ggplot(data = df.10_30, aes(y = TotalSoilC, x = success)) + geom_boxplot(aes(fill = success))

df1 <- read.csv("Data/MTNYCData_modified.csv")

df1$coresection <- 3
i <- 1
df1$coresection <- 3
for (i in 1:392) {
  if (df1$CoreSection_cm[i] == '0-10') {
    df1$coresection[i] = 1
  }
  if (df1$CoreSection_cm[i] == '10--30') {
    df1$coresection[i] = 2
  }
}
write.csv(df1, "Data/MTNYCData_updated.csv")

for (i in 1:392) {
  if (df1$CoreSection_cm[i] == '0-10') {
    df1$coresection[i] = 1
  }
  if (df1$CoreSection_cm[i] == '10--30') {
    df1$coresection[i] = 2
  }
  if (df1$CoreSection_cm[i] == '30-70') {
    df1$coresection[i] = 3
  }
  if (df1$CoreSection_cm[i] == '70-90') {
    df1$coresection[i] = 4
  }
  if (df1$CoreSection_cm[i] == '90-100') {
    df1$coresection[i] = 5
  }
}

library(dplyr)
df_list <- df1 %>% 
  group_split(coresection)

reg_model <- formula(TotalSoilC ~ (1|Site) + factor(Month) + factor(Success) + factor(Diversity) + RootMass_g+Moisture_g+Respiration)

reg_results <- lapply(df_list, function(x) {
  lmer(reg_model, data = x)
})

reg_summary <- lapply(reg_results, summary)

reg_summary[[2]]

## ggpairs to check the multilinearity vif
## different scale: use z score to transform all the numeric independent variables
## + factor(section) + (1|Plot); or separately analyze 3 datasets; or multiple the section to other vars
## y distributed continuously, need to test normality
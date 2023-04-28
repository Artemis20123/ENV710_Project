# to do the subgroup regressions at the same time

library(dplyr)
library(lme4)
library(lmerTest)
library(lmtest)

df1 <- read.csv("Data/MTNYCData_modified2_noDEAna.csv")
# df1 <- read.csv("Data/MTNYCData_modified.csv")

# dept vars normality
shapiro.test(df1$BiomassC)
shapiro.test(log(df1$BiomassC))

# model reduction
lm1 <- lmer(log(BiomassC) ~ Success + factor(diversity) + factor(Season) + Respiration + BiomassN + NO2_NO3 + NH4 + RootMass_g + Moisture_g + factor(coresection) + Year +
                (1|Site/Plot), data = df1)
summary(lm1)

lm2 <- update(lm1, ~.-factor(Season))
summary(lm2)

lm3 <- update(lm2, ~.-NH4)
summary(lm3)

lm4 <- update(lm3, ~.-RootMass_g)
summary(lm4)

lm5 <- update(lm4, ~.-Year)
summary(lm5)

lm6 <- update(lm5, ~.-factor(diversity))
summary(lm6)

AIC(lm1, lm2, lm3,lm4,lm5,lm6)

# assumption checks
plot(lm6)

mean(residuals(lm6))
shapiro.test(residuals(lm6))

vif(lm6)

pchisq(lmm6$deviance, lmm6$df.residual, lower=F)

library(MuMIn)
MuMIn::r.squaredGLMM(lm6)

library(merTools)
merTools::plotFEsim(FEsim(lm6), intercept = T)
library(lattice)
dotplot(ranef(lm6))


df_list <- df1 %>% 
  group_split(coresection)

# reg_model <- formula(BiomassC ~ (1|Site) + factor(Month) + factor(Success) + factor(Diversity) + RootMass_g+Moisture_g+Respiration)
reg_model <- formula(BiomassC ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + NH4 + Nitrification + RootMass_g + Moisture_g + Year +
                       (1|Site/Plot))

noNA1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + Year +
                (1|Site/Plot), data = noNAdata)

reg_results <- lapply(df_list, function(x) {
  lmer(reg_model, data = x)
})

reg_summary <- lapply(reg_results, summary)

reg_summary[[1]]

success.s1.glmer1 <- lme4::glmer(success~(1|Site/Plot)+BiomassN+RockMass_g+Moisture_g+BulkDensity_g_cm3, family=binomial, data=tree_0_10)
summary(success.s1.glmer1)

hist(df1$BiomassN)
shapiro.test(log(df1$BiomassN+7))

library(GGally)
ggpairs(df1, c(3,4,7,8,10,11,14,16))

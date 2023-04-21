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
library(MuMIn)
library(lme4)

#data exploration 

#We need to create a dataset with everything but create a new column that makes soil depth into 
#a factor (low, medium, high). 

#Check correlation and remove correlated variables 

#Still may need to standardize scales. 
#One interaction: soil depth with some other variable. 

modified1 <- read.csv("MTNYCData_modified1.csv")

#screen for correlation 
#source("screen_cor.R")
colnames(modified1)
screen.cor(modified1[,-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 18, 19, 20, 24, 26, 29, 30, 31, 32, 35)])


#Stepwise reduction for BiomassN

m1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + 
              NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + Year + 
              (1|Site/Plot), data = modified1) 
#+ (1|Year)
summary(m1)
plot(m1)
vif(m1)
isSingular(m1)

m2 <- update(m1,.~.-Success)
  
#+ (1|Year)
summary(m2)

#singularity tests 
t3 <- update(m2,.~.-Year)
summary(t3)
t4 <- update(t3,.~.-factor(diversity))
summary(t4)
t5 <- update(t4,.~.-factor(coresection))
summary(t5)

AIC(m1, m2, t3, t4, t5)

m3 <- update(m2,.~.-RootMass_g)
summary(m3)

m4 <- update(m3,.~.-factor(coresection))
summary(m4)
isSingular(m4)
vif(m4)

AIC(m1, m2, m3, m4)

m5 <- update(m4,.~.-factor(diversity))
summary(m5)

isSingular(m5)
AIC(m1, m2, m3, m4, m5)

m6 <- update(m5,.~.-Year)
summary(m6)
AIC(m1, m2, m3, m4, m5, m6)


#TESTING 
TEST <- lmer(BiomassN ~ factor(Season) + Respiration + NO2_NO3 + 
               NH4 + Nitrification + Moisture_g + factor(coresection) + 
             (1|Site) + (1|Plot) + (1|Year), data = modified1)

#finding what variables are causing singularity. 

# calculate the matrix rank and determinant
#rankm4 <- qr(m4)$rank
#detm4 <- det(m4)

# print the rank and determinant
#print(paste0("Matrix rank: ", rankX))
#print(paste0("Matrix determinant: ", detX))



#m4 is all significant 
plot(m4)

#Histograms with the noNA data 

ggplot(aes(x = RootMass_g), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(RootMass_g)), data = noNAdata) + 
  geom_histogram() #NaNs Produced 

ggplot(aes(x = RockMass_g), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(RockMass_g)), data = noNAdata) + 
  geom_histogram() #NaNs Produced 

ggplot(aes(x = Respiration), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(Respiration)), data = noNAdata) + 
  geom_histogram() #Works 

ggplot(aes(x = BiomassN), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(BiomassN)), data = noNAdata) + 
  geom_histogram() #NaNs Produced

ggplot(aes(x = Moisture_g), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(Moisture_g)), data = noNAdata) + 
  geom_histogram() #works 

ggplot(aes(x = BulkDensity_g_cm3), data = noNAdata) + 
  geom_histogram() #fine without log 

ggplot(aes(x = BiomassC), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(BiomassC)), data = noNAdata) + 
  geom_histogram()#Works 

ggplot(aes(x = NO2_NO3), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(NO2_NO3)), data = noNAdata) + 
  geom_histogram() #NaNs produced 

ggplot(aes(x = NH4), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(NH4)), data = noNAdata) + 
  geom_histogram()#NaNs Produced  

ggplot(aes(x = Nitrification), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(Nitrification)), data = noNAdata) + 
  geom_histogram() #NaNs produced 

ggplot(aes(x = DEA), data = noNAdata) + 
  geom_histogram()
ggplot(aes(x = log(DEA)), data = noNAdata) + 
  geom_histogram() #NaNs Produces 

#Start without rootmass at all 

nr1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + 
             NH4 + Nitrification + Moisture_g + factor(coresection) + Year + 
             (1|Site/Plot), data = modified1) 

#Updated dataset with no NA 

noNAdata <- read.csv("MTNYCData_modified2_noDEAna.csv")

noNA1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + 
             NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + Year + 
             (1|Site/Plot), data = noNAdata) 
summary(noNA1)

noNA2 <- update(noNA1,.~.-Success)
summary(noNA2)

noNA3 <- update(noNA2,.~.-Year)
summary(noNA3)

noNA4 <- update(noNA3,.~.-factor(diversity))
summary(noNA4)

noNA5 <- update(noNA4,.~.-RootMass_g)
summary(noNA5)

#Updated dataset with no NA AND log transformed variables dependent variable.
#only worked when only respriration was logged. logging Moisture greated singluarity 

noNAdata <- read.csv("MTNYCData_modified2_noDEAna.csv")
#Respirations, Moisture_g, 
logNoNA1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + log(Respiration) + NO2_NO3 + 
                NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + Year + 
                (1|Site/Plot), data = noNAdata) 
summary(noNA1)

logNoNA2 <- update(logNoNA1,.~.-Success)
summary(logNoNA2)

logNoNA3 <- update(logNoNA2,.~.-Year)
summary(logNoNA3)

logNoNA4 <- update(logNoNA3,.~.-log(Respiration))
summary(logNoNA4)

logNoNA5 <- update(logNoNA4,.~.-RootMass_g)
summary(logNoNA5)

logNoNA6 <- update(logNoNA5,.~.-factor(diversity))
summary(logNoNA6)

#test of all models 
AIC(noNA1, noNA2, noNA3, noNA4, noNA5, logNoNA1, logNoNA2, logNoNA3, logNoNA4, logNoNA5, logNoNA6)

tree <- read.csv("MTNYCData_modified2_noDEAna.csv")
tree.l <- tree
for (i in 1:326){
  tree.l$RockMass_g[i] <- log(tree$RockMass_g[i] + abs(min(tree$RockMass_g)) + 0.001)
  tree.l$RootMass_g[i] <- log(tree$RootMass_g[i] + abs(min(tree$RootMass_g)) + 0.001)
  tree.l$Moisture_g[i] <- log(tree$Moisture_g[i] + abs(min(tree$Moisture_g)) + 0.001)
  tree.l$BulkDensity_g_cm3[i] <- log(tree$BulkDensity_g_cm3[i] + abs(min(tree$BulkDensity_g_cm3)) + 0.001)
  tree.l$BiomassC[i] <- log(tree$BiomassC[i] + abs(min(tree$BiomassC)) + 0.001)
  tree.l$Respiration[i] <- log(tree$Respiration[i] + abs(min(tree$Respiration)) + 0.001)
  tree.l$NO2_NO3[i] <- log(tree$NO2_NO3[i] + abs(min(tree$NO2_NO3)) + 0.001)
  tree.l$NH4[i] <- log(tree$NH4[i] + abs(min(tree$NH4)) + 0.001)
  tree.l$TIN[i] <- log(tree$TIN[i] + abs(min(tree$TIN)) + 0.001)
  tree.l$BiomassN[i] <- log(tree$BiomassN[i] + abs(min(tree$BiomassN)) + 0.001)
  tree.l$DEA[i] <- log(tree$DEA[i] + abs(min(tree$DEA)) + 0.001)
}

no.out.tree.l <- tree.l[-c(33,29,193),]

ggpairs(tree.l[,c(9:17)])
ggpairs(tree.l[,c(18:25)])

#Model construction with logged variables 

finalm1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + 
                   NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + Year + BiomassC +
                   (1|Site/Plot), data = tree.l) 
summary(finalm1)

finalm2 <- update(finalm1,.~.-Year)
summary(finalm2)

finalm3 <- update(finalm2,.~.-Respiration)
summary(finalm3)

finalm4 <- update(finalm3,.~.-RootMass_g)
summary(finalm4)

finalm5 <- update(finalm4,.~.-Moisture_g)
summary(finalm5)

finaFinalModel <- update(finalm5,.~.-Success)
summary(finaFinalModel)

finalm7 <- update(finaFinalModel,.~.-factor(coresection))
summary(finalm7)

#Final model: factor(diversity), factor(season), NO2_NO3, NH4, Nitrification, Biomass C

AIC(finalm1, finalm2, finalm3, finalm4, finalm5, finaFinalModel, finalm7)
plot(finalm7)
vif(finalm7)
check_autocorrelation(finalm7)

colnames(tree.l)
screen.cor(tree.l[,-c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 17, 19, 21, 22, 23, 14, 25)])

#no outliers

out.finalm1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + 
                  NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + Year + BiomassC +
                  (1|Site/Plot), data = no.out.tree.l) 
summary(out.finalm1)

out.finalm2 <- update(out.finalm1,.~.-Respiration)
summary(out.finalm2)

out.finalm3 <- update(out.finalm2,.~.-factor(Season))
summary(out.finalm3)

out.finalm4 <- update(out.finalm3,.~.-Moisture_g)
summary(out.finalm4)

out.finalm5 <- update(out.finalm4,.~.-factor(coresection))
summary(out.finalm5)

out.finaFinalModel <- update(out.finalm5,.~.-Success)
summary(out.finaFinalModel)

out.finalm7 <- update(out.finaFinalModel,.~.-RootMass_g)
summary(out.finalm7)

out.finalm8 <- update(out.finalm7,.~.-Year)
summary(out.finalm8)

AIC(out.finalm1, out.finalm2, out.finalm3, out.finalm4, out.finalm5, out.finaFinalModel, out.finalm7, out.finalm8)
plot(out.finalm7)
vif(out.finalm7)
check_autocorrelation(out.finalm8)

car::influencePlot(out.finalm7)

#Year as a random effect. 

noY.finalm1 <- lmer(BiomassN ~ Success + factor(diversity) + factor(Season) + Respiration + NO2_NO3 + 
                      NH4 + Nitrification + RootMass_g + Moisture_g + factor(coresection) + BiomassC +
                      (1|Site/Plot), data = no.out.tree.l) 
summary(noY.finalm1)

noY.finalm2 <- update(noY.finalm1,.~.-Respiration)
summary(noY.finalm2)

noY.finalm3 <- update(noY.finalm2,.~.-Success)
summary(noY.finalm3)

noY.finalm4 <- update(noY.finalm3,.~.-factor(coresection))
summary(noY.finalm4)

noY.finalm5 <- update(noY.finalm4,.~.-Moisture_g)
summary(noY.finalm5)

noY.finaFinalModel <- update(noY.finalm5,.~.-RootMass_g)
summary(noY.finaFinalModel)

AIC(out.finalm1, out.finalm2, out.finalm3, out.finalm4, out.finalm5, out.finaFinalModel, out.finalm7, out.finalm8,
    noY.finalm1, noY.finalm2, noY.finalm3, noY.finalm4, noY.finalm5, noY.finaFinalModel)
plot(noY.finaFinalModel)
vif(noY.finaFinalModel)
check_autocorrelation(noY.finaFinalModel)

#noY.finalm7 <- update(noY.finaFinalModel,.~.-RootMass_g)
#summary(noY.finalm7)

#FINAL MODEL 

FinalModel <- update(noY.finalm5,.~.-RootMass_g)
summary(FinalModel)
plot(FinalModel)
var.test(FinalModel)
vif(FinalModel)
check_autocorrelation(FinalModel)
r.squaredGLMM(FinalModel)

success.Final.Model <- simulateResiduals(fittedModel = FinalModel)
plot(success.Final.Model)

#BiomassC Plot 
 
Ncoef <- fixef(FinalModel)
Ncoef

success.eq1 <- function(x){Ncoef[1]+Ncoef[2]+
    Ncoef[3]+Ncoef[4]*mean(no.out.tree.l$NO2_NO3)+
    Ncoef[5]*mean(no.out.tree.l$NH4)+Ncoef[6]*mean(no.out.tree.l$Nitrification)+Ncoef[7]*x}
success.eq2 <- function(x){Ncoef[1]+
    Ncoef[3]+Ncoef[4]*mean(no.out.tree.l$NO2_NO3)+
    Ncoef[5]*mean(no.out.tree.l$NH4)+Ncoef[6]*mean(no.out.tree.l$Nitrification)+Ncoef[7]*x}

ggplot(no.out.tree.l, aes(x=BiomassC, y=BiomassN)) + 
  geom_jitter(height = 0.02, width = 0.02, alpha = 0.3, size = 3) +
  stat_function(fun=success.eq1, geom="line", linewidth=0.7, aes(color="High")) +
  stat_function(fun=success.eq2, geom="line", linewidth=0.7, aes(color="Low")) + 
  theme_bw() +
  theme(legend.position = "right") +
  labs(y = "BiomassN", x = "BiomassC") +
  scale_color_viridis("Diversity", discrete = TRUE)



## Assumption checks

### dept var normality


# dept vars normality
shapiro.test(no.out.tree.l$BiomassN)
shapiro.test(log(no.out.tree.l$BiomassN))


### Residuals mean & normality

mean(residuals(FinalModel))
shapiro.test(residuals(FinalModel))

# i think the one below is better
library(DHARMa)
lmm6 <- simulateResiduals(fittedModel = FinalModel) 
plot(lmm6)

library(performance)
performance::check_model(FinalModel, check = c("reqq", "qq", "linearity"))

### Goodness of fit / R2 FE RE

library(MuMIn)
MuMIn::r.squaredGLMM(FinalModel)

library(merTools)
merTools::plotFEsim(FEsim(FinalModel), intercept = T)
library(lattice)
dotplot(ranef(FinalModel))

#In the random effect term we are indicating that campaign is nested within city.

### Outlier test

library(car)
outlierTest(FinalModel)
influencePlot(FinalModel)

### Multicollinearity

vif(FinalModel)
performance::check_model(FinalModel, check = c("vif"))




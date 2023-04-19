# to do the subgroup regressions at the same time

library(dplyr)

df1 <- read.csv("Data/MTNYCData_modified.csv")

df_list <- df1 %>% 
  group_split(coresection)

reg_model <- formula(TotalSoilC ~ (1|Site) + factor(Month) + factor(Success) + factor(Diversity) + RootMass_g+Moisture_g+Respiration)

reg_results <- lapply(df_list, function(x) {
  lmer(reg_model, data = x)
})

reg_summary <- lapply(reg_results, summary)

reg_summary[[1]]

library(GGally)
ggpairs(df1)

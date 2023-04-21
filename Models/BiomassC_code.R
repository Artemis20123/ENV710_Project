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
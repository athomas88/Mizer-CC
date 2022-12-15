# Calculate upper and lower bounds for sensitivity analysis -------------------------------------------------------------------------
# description
library(tidyverse)
library(fishmethods)
setwd("C:/Users/athomas/Documents/Mizer_cc/")

mizer_spp <- read_csv("data/myanmar_model_inputs_54_TL.csv",col_names=TRUE)
spp_db <- read_csv("C:/Users/athomas/Documents/Myanmar/Mizer/data/Species_FB_SLB_M_TL_v2.csv",col_names = TRUE)
grps <- unique(mizer_spp$FamilyGroup)

grpBounds <- data.frame("FamilyGroup"=character(),
                        "Linf_low"=double(),
                        "Linf_high"=double(),
                        "K_low"=double(),
                        "K_high"=double(),
                        "Lmat_low"=double(),
                        "Lmat_high"=double(),
                        "TL_low"=double(),
                        "TL_high"=double())

for (i in grps){
  grp <- filter(spp_db, FamilyGroup == i)
  bounds <- data.frame(i,
                       min(grp$L_inf_cm, na.rm = TRUE),
                       max(grp$L_inf_cm, na.rm = TRUE),
                       min(grp$VB_k, na.rm = TRUE),
                       max(grp$VB_k, na.rm = TRUE),
                       min(grp$L_maturity, na.rm = TRUE),
                       max(grp$L_maturity, na.rm = TRUE),
                       min(grp$trophic_level, na.rm = TRUE),
                       max(grp$trophic_level, na.rm = TRUE))
  names(bounds) <- names(grpBounds)
  grpBounds <- rbind(grpBounds,bounds)
}

spp_sens <- merge(mizer_spp, grpBounds, by="FamilyGroup")
spp_sens[spp_sens=="Inf"] <- NA
spp_sens[spp_sens=="-Inf"] <- NA

#Calculate new Ms
#spp_sens <- read_csv("data/myanmar_model_inputs_54_sensitivity.csv",col_names=TRUE)

for (i in 1:nrow(spp_sens)){
  spp_sens$M_lowL[i] <- as.numeric(M.empirical(Linf = spp_sens$Linf_low[i], Winf = NULL, Kl = spp_sens$VB_k[i], Kw = NULL,
                                               TC = spp_sens$Temperature[i], tmax = NULL, tm = NULL, GSI = NULL, Wdry = NULL,
                                               Wwet = NULL, Bl = NULL, method = c(11)))
  spp_sens$M_highL[i] <- as.numeric(M.empirical(Linf = spp_sens$Linf_high[i], Winf = NULL, Kl = spp_sens$VB_k[i], Kw = NULL,
                                               TC = spp_sens$Temperature[i], tmax = NULL, tm = NULL, GSI = NULL, Wdry = NULL,
                                               Wwet = NULL, Bl = NULL, method = c(11)))
  spp_sens$M_lowK[i] <- as.numeric(M.empirical(Linf = spp_sens$L_inf_cm[i], Winf = NULL, Kl = spp_sens$K_low[i], Kw = NULL,
                                                TC = spp_sens$Temperature[i], tmax = NULL, tm = NULL, GSI = NULL, Wdry = NULL,
                                                Wwet = NULL, Bl = NULL, method = c(11)))
  spp_sens$M_highK[i] <- as.numeric(M.empirical(Linf = spp_sens$L_inf_cm[i], Winf = NULL, Kl = spp_sens$K_high[i], Kw = NULL,
                                               TC = spp_sens$Temperature[i], tmax = NULL, tm = NULL, GSI = NULL, Wdry = NULL,
                                               Wwet = NULL, Bl = NULL, method = c(11)))

}
write_csv(spp_sens,"data/myanmar_model_inputs_54_sensitivity.csv")

#Manually replaced NA values with the spp value after writing
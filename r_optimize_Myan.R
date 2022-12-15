# global.r
## Myanmar multispecies finfish study
## Calculate swept area biomass and set up optimization of biomass to calculate r.
## Kristin Kleisner
## June 15, 2018
## Load packages
rm(list = ls())
library(mizer); library(tidyverse); library(rhandsontable); library(gridExtra)
#library(rfishbase); library(fishmethods); library(magrittr); library(ggplot2)

setwd('C:/Users/kkleisner/Documents/ContractWork/EDF/Upside/Helmsley/Myanmar/')
 N_Bio<-read.csv('Nansen_Biomass_SweptArea.csv', header=T, stringsAsFactors = F)
 # Check<-read.csv('Check_Genus_Group_Calib.csv', header=T, stringsAsFactors = F)
 # 
 # Check2<-Check%>%
 #   group_by(Genus)%>%
 #   summarize(aveLinf=mean(L_inf_cm, na.rm=T), sdLinf=sd(L_inf_cm, na.rm = T), ctLinf=length(L_inf_cm), 
 #             avetmax=mean(tmax, na.rm=T), sdtmax=sd(tmax, na.rm = T),
 #             aveLmat=mean(L_maturity, na.rm=T), sdLmat=sd(L_maturity, na.rm = T),
 #             avetm=mean(tm, na.rm=T), sdtm=sd(tm, na.rm = T))
 
##Calculate swept area: 
#SArea = velocity (V)*time (t)*headrope length (hr)*width of path swept by trawl (X2)
#This equation is from Sparre and Venema (1998), an FAO report
#According to the 2013 Nansen report, hr*X2: 
#"The effective fishing width of trawl gear used by R/V Dr Fridtjof Nansen is considered to be 18.5 m."
#Need to convert duration from minutes to hours. Velocity is given in knots (nm/hr):
N_Bio<-N_Bio%>%
  mutate(Dur_hr=Duration/60, SArea=Speed*Dur_hr*18.5, CW_a=CH_weight/(SArea/Dur_hr))

N_Bio2<-N_Bio%>%
  group_by(Region)%>%
  summarize(TotalArea=sum(SArea, na.rm=T))
N_Bio3<-N_Bio%>%
  full_join(N_Bio2)%>%
  group_by(Region, sciName, Mizer, TotalArea)%>%
  summarize(CW_a_spp=sum(CW_a, na.rm=T))%>%
  mutate(current_biomass_kg=CW_a_spp*TotalArea)

##There is no Lutjanus lutjanus biomass in Rakhine, so we will use the Tanintharyi estimate:
N_Bio3_LL<-N_Bio3[which(N_Bio3$sciName=="Lutjanus lutjanus"),]


setwd('C:/Users/kkleisner/Documents/ContractWork/EDF/Upside/Helmsley/Myanmar/myanmar-mizer/')
source("functions.R")
## Load data:
params_site_specific <- read_csv("data/params_myanmar_site_specific_TL.csv",col_names=TRUE)
params_site_specific <- params_site_specific %>%
  left_join(N_Bio3)
params_site_specific$current_biomass_kg<-ifelse(params_site_specific$sciName=="Lutjanus lutjanus", N_Bio3_LL$current_biomass_kg, params_site_specific$current_biomass_kg)

## Massage data
params_default <- params_site_specific %>%
  filter(!is.na(commonName)) %>%
  mutate(w_a = unitConversion(unit,Weight_a_cm_g,Weight_b_cm_g,FL_TL_a,FL_TL_b,SL_TL_a,SL_TL_b),
         w_b = Weight_b_cm_g,
         w_inf = w_a * L_inf_cm ^ w_b,
         w_mat = w_a * L_maturity_cm ^ w_b,
         beta = 100,
         sigma = 2,
         r_max = 1e9 * Productivity /3,
         k_vb = VB_k,
         species = Stock,
         M = natural_mortality_M,
         FvM = FvM_current)  %>%
  dplyr::select(species,w_inf,w_mat,beta,sigma,r_max,k_vb,w_a,w_b,M,Productivity,trophic_level,FvM,current_biomass_kg) %>%
  as.data.frame()
## Set up mizer parameters
## Define interactions matrix. Assume all species spatially overlap and interact with each other
## So interactions are based purely on size
inter_spp <- matrix(1,nrow=nrow(params_default),ncol=nrow(params_default))
rownames(inter_spp) <- params_default$species
colnames(inter_spp) <- params_default$species
## set up community and background spectrum
#tVirgin <- 1800 # year to start virgin simulation
#t0 <- 1900 # year to start historic simulation
tManagement <- 2017 # year to start forward projection, when new management starts
#tEnd <- 2067 # year to end forward projection
nsS <- nrow(params_default) ## number of species
nwS <- 100 ## number of size bins in community spectrum
npS <- 30 ## extra size bins in background spectrum
kappa <- 9e12 ## resource spectrum



optimizationFunction <- function(x){
params_default <- params_default %>%
  mutate(r_max = x * Productivity / 3)

MizerParams <- MizerParams(params_default,no_w=nwS,no_w_pp=npS,interaction = inter_spp,kappa = kappa)
MizerParams@species_params$M <- params_default$M

tRes <- 1
t0 <-tManagement - 50
tVirgin <- t0 - 50
tEnd <- tManagement + 25
fisheryAge <- 50
## Define virgin unfished period
virginInput <- list(managementName = "Virgin Unfished",
                    interactions = inter_spp,
                    params_data = params_default,
                    sel_func = "knife_edge",
                    knife_edge_size = 0,
                    time = seq(tVirgin,t0-1),
                    e0 = rep(0,length(MizerParams@species_params$M)),
                    eDelta = 0,
                    initN = get_initial_n(MizerParams),
                    initNPP = MizerParams@cc_pp,
                    nws = nwS,
                    nps = npS,
                    nss = nsS,
                    res = tRes)
virginOutput <- tidySim(virginInput)
virginN <- virginOutput$sim@n[dim(virginOutput$sim@n)[1],,]
virginNPP <- virginOutput$sim@n_pp[dim(virginOutput$sim@n_pp)[1],]
#progress$set(value = progressCounter, detail = paste(round(progressCounter*100/progressTotal,0),"% done.",sep=""))
## Define historic fishing period
FvM = params_default$FvM
historicInput <- list(managementName = "Historic Fishing",
                      interactions = inter_spp,
                      params_data = params_default,
                      sel_func = "knife_edge",
                      knife_edge_size = 0,
                      time = seq(t0,tManagement-1),
                      e0 = FvM/fisheryAge,
                      eDelta = (fisheryAge)^(1/fisheryAge)-1,
                      initN = virginN,
                      initNPP = virginNPP,
                      nws = nwS,
                      nps = npS,
                      nss = nsS,
                      res = tRes)
historicOutput <- tidySim(historicInput)

biomassIndicator <- historicOutput$simResultsAggregated %>%
  filter(Year == tManagement -1) %>%
  select(Biomass) %>%
  as.numeric() /1000

biomassReference <- sum(params_default$current_biomass_kg)

diff <- abs(biomassIndicator - biomassReference)
print(paste("r_max guess: ",x,sep=""))
print(paste("Difference: ",diff,sep=""))
if(diff<1) diff <- 0
return(diff)
}

## ADD code to create starting guess, average of the r_max of all unique spp in model

#starting guess (average of 1e9*Productivity/3 for 22 species): 8.416666667e8
OptR <- optim(par=8.416666667e8,optimizationFunction,upper=3e11,lower=1e7,method="L-BFGS-B",control=list(maxit=1000))
save(OptR,file="data/MyanOptR.Rdata") #"Difference: 68362233.6303182"

## ADD code to calculate new r_max for each spp based on the OptR value : OptR value * Productivity / 3 or OptR value * Trophic Level / 4

#load('data/MyanOptR.Rdata')

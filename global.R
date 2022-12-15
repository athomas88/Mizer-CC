# global.r
## Myanmar multispecies finfish study
## October 18, 2017
##Original code created by Gavin McDonald and modified by Kristin Kleisner
## Load packages
rm(list = ls())
library(mizer)
library(tidyverse)
library(rhandsontable)
library(gridExtra)
library(rstudioapi)
library(shiny)
library(shinyWidgets)

#Set the directory to the current location to access scripts and data
source("functions.R")
source("CC_functions.R")
## Load data
params_belize_default <- read_csv("data/myanmar_model_inputs_gm.csv",col_names=TRUE)
#params_belize_site_specific <- read_csv("data/params_belize_site_specific.csv",col_names=TRUE)
## Massage data
params_belize_default <- params_belize_default %>%
  filter(!is.na(Stock)) %>%
  mutate(w_a = unitConversion(unit,Weight_a_cm_g,Weight_b_cm_g,FL_TL_a,FL_TL_b,SL_TL_a,SL_TL_b),
         w_b = Weight_b_cm_g,
         w_inf = w_a * L_inf_cm ^ w_b,
         w_mat = w_a * L_maturity_cm ^ w_b,
         beta = 100,
         sigma = 2,
         r_max = r_max,
         k_vb = VB_k,
         sciname = sciName,
         species = Stock,
         region = Region,
         M = natural_mortality_M,
         FvM = F  )%>%
  dplyr::select(species, sciname, commonName, region, w_inf,w_mat,beta,sigma,k_vb,w_a,w_b,M,trophic_level,Vulnerability, 
                Vuln_Cheung, Resilience, PriceCat,FvM,r_max) %>%
  as.data.frame()
## Set up mizer parameters
## Define interactions matrix. Assume all species spatially overlap and interact with each other
## So interactions are based purely on size
inter_belize <- matrix(1,nrow=nrow(params_belize_default),ncol=nrow(params_belize_default))
rownames(inter_belize) <- params_belize_default$species
colnames(inter_belize) <- params_belize_default$species
## set up community and background spectrum
#tVirgin <- 1800 # year to start virgin simulation
#t0 <- 1900 # year to start historic simulation
tNow <- 2019 # year to start forward projection, when new management starts
#tEnd <- 2067 # year to end forward projection
nsS <- nrow(params_belize_default) ## number of species
nwS <- 100 ## number of size bins in community spectrum
npS <- 30 ## extra size bins in background spectrum
kappa <- 9e12 ## resource spectrum
fisheryAge <- 50
t0 <-tNow - fisheryAge

sub_rakhine<-params_belize_default%>%
  filter(region=="Rakhine")
sub_delta<-params_belize_default%>%
  filter(region=="Delta")
sub_tanin<-params_belize_default%>%
  filter(region=="Tanintharyi")

monthlyCatch <- read.csv("data/SeasonalCatch.csv")
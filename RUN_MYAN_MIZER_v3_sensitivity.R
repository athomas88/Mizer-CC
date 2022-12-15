# global.r
## Myanmar multispecies finfish study
## Load packages

rm(list = ls())

setwd("C:/Users/athomas/Documents/Mizer_cc/")
## Load data

library(mizer)
library(tidyverse)
source("functions.R")
library(rhandsontable)
library(gridExtra)


params_belize_default <- read_csv("data/myanmar_model_inputs_54_sensitivity.csv",col_names=TRUE)
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
         Stock = Stock,
         species = sciName,
         region = Region,
         M = natural_mortality_M,
         tropic_level = trophic_level,
         FvM = F1)%>%
  dplyr::select(species, Stock, commonName, region, w_inf,w_mat,beta,sigma,k_vb,w_a,w_b,M,trophic_level,Vulnerability, 
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
tManagement <- 2014 # year to start forward projection, when new management starts
#tEnd <- 2067 # year to end forward projection
nsS <- nrow(params_belize_default) ## number of species
nwS <- 100 ## number of size bins in community spectrum
npS <- 30 ## extra size bins in background spectrum
kappa <- 9e12 ## resource spectrum

sub_rakhine<-params_belize_default%>%
  filter(region=="Rakhine")
sub_delta<-params_belize_default%>%
  filter(region=="Delta")
sub_tanin<-params_belize_default%>%
  filter(region=="Tanintharyi")
sub_myanmar<-params_belize_default%>%
  filter(region=="Myanmar")



#################################################################

#values$params_belize<-values$params_temp
#progress <- shiny::Progress$new()
#progress$set(message = "Running simulation:", value = 0.01)
#on.exit(progress$close())
#progressCounter <- 0
#progressTotal <- 2 + length(input$interventionGroup)
# if (input$siteAnalysis == "All"){
#   params_belize_analysis <- values$params_belize
# }
# if (input$siteAnalysis == "Rakhine"){
#   params_belize_analysis <- values$params_belize %>%
#     filter(region == "Rakhine")
# }
# if (input$siteAnalysis == "Delta"){
#   params_belize_analysis <- values$params_belize %>%
#     filter(region == "Delta")
# }
# if (input$siteAnalysis == "Tanintharyi"){
#   params_belize_analysis <- values$params_belize %>%
#     filter(region == "Tanintharyi")
# }
params_belize_analysis <-sub_myanmar
values<-list()

nsS <- nrow(params_belize_analysis)
inter_belize <- matrix(1,nrow=nrow(params_belize_analysis),ncol=nrow(params_belize_analysis))
rownames(inter_belize) <- params_belize_analysis$species
colnames(inter_belize) <- params_belize_analysis$species
values$params_belize<-params_belize_analysis

MizerParamsBelize <- MizerParams(params_belize_analysis,no_w=nwS,no_w_pp=npS,interaction = inter_belize,kappa = kappa)
MizerParamsBelize@species_params$M <- params_belize_analysis$M

###THIS input TAKES THE PLACE OF THE UI INPUTS:
input<-NULL
input$timeStep<-0.1
input$fisheryAge<-50
input$tProjection<-86
input$effortCreep<-1
input$singleF<-0.5
input$ntzSize<-20
input$singleSizeLimit<-10

weightLimits <- params_belize_analysis$w_a * input$singleSizeLimit ^ params_belize_analysis$w_b

tRes <- input$timeStep
t0 <-tManagement - input$fisheryAge
tVirgin <- t0 - 50
tEnd <- tManagement + input$tProjection
## Define virgin unfished period
virginInput <- list(managementName = "Virgin unfished",
                    interactions = inter_belize,
                    params_data = params_belize_analysis,
                    sel_func = "knife_edge",
                    knife_edge_size = 0,
                    time = seq(tVirgin,t0-1),
                    e0 = rep(0,length(MizerParamsBelize@species_params$M)),
                    eDelta = 0,
                    initN = get_initial_n(MizerParamsBelize),
                    initNPP = MizerParamsBelize@cc_pp,
                    nws = nwS,
                    nps = npS,
                    nss = nsS,
                    res = tRes)
#progress$set(detail = virginInput$managementName)
virginOutput <- tidySim(virginInput)
virginN <- virginOutput$sim@n[dim(virginOutput$sim@n)[1],,]
virginNPP <- virginOutput$sim@n_pp[dim(virginOutput$sim@n_pp)[1],]
#progressCounter <- progressCounter + 1
#progress$set(value = progressCounter, detail = paste(round(progressCounter*100/progressTotal,0),"% done.",sep=""))
## Define historic fishing period
historicInput <- list(managementName = "Historic fishing",
                      interactions = inter_belize,
                      params_data = params_belize_analysis,
                      sel_func = "knife_edge",
                      knife_edge_size = 0,
                      time = seq(t0,tManagement-1),
                      e0 = params_belize_analysis$FvM/input$fisheryAge,
                      eDelta = (input$fisheryAge)^(1/input$fisheryAge)-1,
                      initN = virginN,
                      initNPP = virginNPP,
                      nws = nwS,
                      nps = npS,
                      nss = nsS,
                      res = tRes)
#progress$set(value = progressCounter/progressTotal,detail = historicInput$managementName)
historicOutput <- tidySim(historicInput)
baselineN <- historicOutput$sim@n[dim(historicOutput$sim@n)[1],,]
baselineNPP <- historicOutput$sim@n_pp[dim(historicOutput$sim@n_pp)[1],]
#progressCounter <- progressCounter + 1
#progress$set(value = progressCounter, detail = paste(round(progressCounter*100/progressTotal,0),"% done.",sep=""))
# if (7 %in% input$interventionGroup){
#   weightLimits <- params_belize_analysis$w_a * input$singleSizeLimit ^ params_belize_analysis$w_b
# } else weightLimits <- vector()
# # if (8 %in% input$interventionGroup){
#   eBasket <- vulnerability(params_belize_analysis$Vulnerability,input$lowBasketF,input$mediumBasketF,input$highBasketF,input$veryHighBasketF)
# } else eBasket <- vector()
# if (9 %in% input$interventionGroup){
#   fTime <- vulnerability(params_belize_analysis$Vulnerability,input$lowBasketF,input$mediumBasketF,input$highBasketF,input$veryHighBasketF)
# } else fTime <- vector()

eBasket <- vector()
eDelta <- input$effortCreep/100
forwardProjectionInputs <- list(
  statusQuo = list(managementName = "No change",
                   interactions = inter_belize,
                   params_data = params_belize_analysis,
                   sel_func = historicInput$sel_func,
                   knife_edge_size = historicInput$knife_edge_size,
                   time = seq(tManagement-1,tEnd),
                   e0 = params_belize_analysis$FvM,
                   eDelta = eDelta,
                   initN = baselineN,
                   initNPP = baselineNPP,
                   nws = nwS,
                   nps = npS,
                   nss = nsS,
                   res = tRes),
  
  FInput = list(managementName = "Limit catch relative to natural mortality",
                interactions = inter_belize,
                params_data = params_belize_analysis,
                sel_func = "knife_edge",
                knife_edge_size = historicInput$knife_edge_size,
                time = seq(tManagement-1,tEnd),
                e0 = params_belize_analysis$M,
                eDelta = 0,
                initN = baselineN,
                initNPP = baselineNPP,
                nws = nwS,
                nps = npS,
                nss = nsS,
                res = tRes),
  
  sizeLimitInput = list(managementName = "Minimum size limit (for each species)",
                        interactions = inter_belize,
                        params_data = params_belize_analysis,
                        sel_func = "knife_edge",
                        knife_edge_size = params_belize_analysis$w_mat,
                        time = seq(tManagement-1,tEnd),
                        e0 = params_belize_analysis$FvM,
                        eDelta = eDelta,
                        initN = baselineN,
                        initNPP = baselineNPP,
                        nws = nwS,
                        nps = npS,
                        nss = nsS,
                        res = tRes),
  
  sizeAndFInput = list(managementName = "Size limits and catch limits",
                       interactions = inter_belize,
                       params_data = params_belize_analysis,
                       sel_func = "knife_edge",
                       knife_edge_size = params_belize_analysis$w_mat,
                       time = seq(tManagement-1,tEnd),
                       e0 = params_belize_analysis$M,
                       eDelta = 0,
                       initN = baselineN,
                       initNPP = baselineNPP,
                       nws = nwS,
                       nps = npS,
                       nss = nsS,
                       res = tRes),
  
  ntz = list(managementName = "No-take zone",
             interactions = inter_belize,
             params_data = params_belize_analysis,
             sel_func = historicInput$sel_func,
             knife_edge_size = historicInput$knife_edge_size,
             time = seq(tManagement-1,tEnd),
             e0 = params_belize_analysis$FvM*(1 - input$ntzSize/100),
             eDelta = eDelta,
             initN = baselineN,
             initNPP = baselineNPP,
             nws = nwS,
             nps = npS,
             nss = nsS,
             res = tRes),
  
  FInputSingle = list(managementName = "Single F across species",
                      interactions = inter_belize,
                      params_data = params_belize_analysis,
                      sel_func = "knife_edge",
                      knife_edge_size = historicInput$knife_edge_size,
                      time = seq(tManagement-1,tEnd),
                      e0 = rep(input$singleF,nsS),
                      eDelta = 0,
                      initN = baselineN,
                      initNPP = baselineNPP,
                      nws = nwS,
                      nps = npS,
                      nss = nsS,
                      res = tRes),
  
  sizeInputSingle = list(managementName = "Single minimum size limit across species",
                         interactions = inter_belize,
                         params_data = params_belize_analysis,
                         sel_func = "knife_edge",
                         knife_edge_size = weightLimits,
                         time = seq(tManagement-1,tEnd),
                         e0 = params_belize_analysis$FvM,
                         eDelta = 0,
                         initN = baselineN,
                         initNPP = baselineNPP,
                         nws = nwS,
                         nps = npS,
                         nss = nsS,
                         res = tRes),
  
  fBasketSingle = list(managementName = "F based on vulnerability basket",
                       interactions = inter_belize,
                       params_data = params_belize_analysis,
                       sel_func = "knife_edge",
                       knife_edge_size = historicInput$knife_edge_size,
                       time = seq(tManagement-1,tEnd),
                       e0 = eBasket,
                       eDelta = 0,
                       initN = baselineN,
                       initNPP = baselineNPP,
                       nws = nwS,
                       nps = npS,
                       nss = nsS,
                       res = tRes))

####ADDED SOME CODE TWEAKS HERE TO MAKE THIS ABLE TO RUN MULTIPLE MANAGEMENT SCENARIOS:  


  forwardProjectionInputs2 <- forwardProjectionInputs[as.numeric(1:7)]
  
  values$forwardProjectionOutputs <- map(forwardProjectionInputs2,function(x){
    #progress$set(value = progressCounter/progressTotal,detail = x$managementName)
    #progressCounter <<- progressCounter + 1
    tidySim(x)
    
  })


values$combinedSpeciesOutput <- values$forwardProjectionOutputs %>%
#values$combinedSpeciesOutput <- FPInputs %>%
  map("simResultsSpecies") %>% 
  bind_rows() %>%
  rbind(historicOutput$simResultsSpecies,
        virginOutput$simResultsSpecies) %>%
  filter(!is.na(Year))

values$combinedAggregatedOutput <- values$forwardProjectionOutputs %>%
  map("simResultsAggregated") %>% 
  bind_rows() %>%
  rbind(historicOutput$simResultsAggregated,
        virginOutput$simResultsAggregated) %>%
  filter(!is.na(Year))

speciesVirgin <- rownames(virginN)
values$virginTidy <- virginN %>% as_data_frame()
values$virginTidy$Species <- speciesVirgin
values$virginTidy <- values$virginTidy %>%
  gather(weightClass,virginAbundance,-Species)
values$virginTidy$weightClass <- as.numeric(values$virginTidy$weightClass)

###################
simV<-virginOutput %>% map("sim")
#virginOutput@n[,"Alectis indica",]
###################


simS<-values$forwardProjectionOutputs %>% map("sim")
timeStep <- input$tProjection +1
simNames <- names(simS)
mapNames <- map(forwardProjectionInputs,function(x){x$managementName}) %>% as.vector()
values$sizeSpectrum <- map(simNames,function(x){
  temp <- get(x,simS)
  mName <- get(x,mapNames)
  temp <- temp@n[timeStep,,]
  speciesName <- rownames(temp)
  temp <- temp %>% as_data_frame()
  temp$Species <- speciesName
  temp <- temp %>%
    gather(weightClass,n,-Species)
  temp$weightClass <- as.numeric(temp$weightClass)
  temp$Management <- mName
  temp <- temp %>%
    mutate(Biomass = weightClass * n) %>%
    filter(n>0) %>%
    left_join(values$virginTidy) %>%
    mutate(relativeAbundance = n/virginAbundance)
}) %>%
  bind_rows()

#### CALCULATE FINAL VALUES
tEnd <- tManagement + input$tProjection
finalValues <- values$combinedAggregatedOutput %>%
  filter(Year == tEnd)

# Merge the outputs from running above code -------------------------------------------------------------------------
# description

finalValues$Parameters <- "Original"
finalValuesLint_low$Parameters <- "Low Linf"
finalValuesLint_high$Parameters <- "High Linf"
finalValuesLmat_low$Parameters <- "Low Lmat"
finalValuesLmat_high$Parameters <- "High Lmat"
finalValuesM_low$Parameters <- "Low M"
finalValuesM_high$Parameters <- "High M"
finalValuesK_low$Parameters <- "Low K"
finalValuesK_high$Parameters <- "High K"


sensitivityValues <- rbind(finalValuesLint_high, finalValuesLint_low)
sensitivityValues <- rbind(sensitivityValues, finalValues)

write_csv(sensitivityValues, "data/sensitivtyOutputs3.csv")

# reformat output for comparison -------------------------------------------------------------------------
# description

#origValues <- read_csv("data/sensitivtyOutputs.csv")

biomass <- sensitivityValues
yield <- sensitivityValues



biomass <- biomass[,c(1,4,5)]
biomass <- pivot_wider(biomass, names_from = Parameters, values_from = Biomass)
write_csv(biomass, "data/sensitivity3_biomass.csv")

yield <- yield[,c(2,4,5)]
yield <- pivot_wider(yield, names_from = Parameters, values_from = Yield)
write_csv(yield, "data/sensitivity3_yield.csv")

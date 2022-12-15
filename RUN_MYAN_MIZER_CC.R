# RUN CC MODEL -------------------------------------------------------------------------
# Myanmar multispecies finfish study
# Currently code is set up to run both non-cc and cc mizers, and then either plot them individually or create comparative plots

# Create workspace/library -------------------------------------------------------------------------
# Necessary setup steps to run code
## Navigate to where you downloaded the GitHub/code
setwd("C:/Users/athomas/Documents/Mizer_cc/")


library(mizer)
library(tidyverse)
source("functions.R")
source("cc_functions.R")
library(rhandsontable)
library(gridExtra)


# Load data -------------------------------------------------------------------------
# Navigate to downloaded spp inputs
params_belize_default <- read_csv("data/myanmar_model_inputs_gm.csv",col_names=TRUE)
#params_belize_site_specific <- read_csv("data/params_belize_site_specific.csv",col_names=TRUE)

# Calculate climate rate -------------------------------------------------------------------------
# Average yearly change in range from 2012 to 2100 by species
ccImpact <- read.csv("data/Myanmar_rangearea_cc.csv")
yrtot <- 88
ccImpact[["cRate"]] <- (ccImpact[["yr00"]] / ccImpact[["yr12"]]) ^ (1 / yrtot)
ccImpact <- ccImpact[,c(1,3,24)]
ccImpact[["Species"]] <- gsub("_", " ", ccImpact[["Species"]])
for (i in 1:nrow(ccImpact)){
  if (ccImpact[["cRate"]][i] >= 1){
    ccImpact[["cRate"]][i] <- 0
  }
  else{
    ccImpact[["cRate"]][i] <- 1 - ccImpact[["cRate"]][i]
  }
}
impact <- ccImpact[ccImpact$RCP==8.5,]
row.names(impact) <- impact[[1]]

# Alter Climate variable -------------------------------------------------------------------------
# Optional to test impact of climate variable
#impact[3] <- 0.5
impact[3] <- impact[3] * 1

# Massage data -------------------------------------------------------------------------
# Transform model inputs into format for running the model
params_belize_default$F1 <- params_belize_default$F
params_belize_default$F1[params_belize_default$F1 > 3] <- 3 #Cap F at 3

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
         FvM = F1  )%>%
  dplyr::select(species, Stock, commonName, region, w_inf,w_mat,beta,sigma,k_vb,w_a,w_b,M,trophic_level,Vulnerability, 
                Vuln_Cheung, Resilience, PriceCat,FvM,r_max) %>%
  as.data.frame()

# Parameter setup -------------------------------------------------------------------------
# Base model parameters

## Define interactions matrix. Assume all species spatially overlap and interact with each other
## So interactions are based purely on size
inter_belize <- matrix(1,nrow=nrow(params_belize_default),ncol=nrow(params_belize_default))
rownames(inter_belize) <- params_belize_default$species
colnames(inter_belize) <- params_belize_default$species

## Set up community and background spectrum
tManagement <- 2014 # year to start forward projection, when new management starts
nsS <- nrow(params_belize_default) ## number of species
nwS <- 100 ## number of size bins in community spectrum
npS <- 30 ## extra size bins in background spectrum
kappa <- 9e12 ## resource spectrum

## Create regional inputs
## spp that are found in multiple regions have multiple rows in the raw data
## (one for each region) so need to run regions separately
## have all-region data for Kristin's 54 spp for later addition
sub_rakhine<-params_belize_default%>%
  filter(region=="Rakhine")
sub_delta<-params_belize_default%>%
  filter(region=="Delta")
sub_tanin<-params_belize_default%>%
  filter(region=="Tanintharyi")


# Set up inputs -------------------------------------------------------------------------
# Get specific info from the region inputs to start model runs
params_belize_analysis <- sub_tanin
values<-list()
valuesCC<-list()

nsS <- nrow(params_belize_analysis)
inter_belize <- matrix(1,nrow=nrow(params_belize_analysis),ncol=nrow(params_belize_analysis))
rownames(inter_belize) <- params_belize_analysis$species
colnames(inter_belize) <- params_belize_analysis$species
values$params_belize<-params_belize_analysis
valuesCC$params_belize<-params_belize_analysis

# Create regional MizerParams -------------------------------------------------------------------------
# This will be the base parameters for all model runs
MizerParamsBelize <- MizerParams(params_belize_analysis,no_w=nwS,no_w_pp=npS,interaction = inter_belize,kappa = kappa)
MizerParamsBelize@species_params$M <- params_belize_analysis$M

# Replace UI inputs -------------------------------------------------------------------------
# These would be defined by the Shiny app, need to manually define them
input<-NULL
input$timeStep<-0.1
input$fisheryAge<-50
input$tProjection<-86
input$effortCreep<-1
input$singleF<-0.5
input$ntzSize<-20
input$singleSizeLimit<-10

# Final parameters -------------------------------------------------------------------------
# Parameters that were dependent on the MizerParams or the user inputs
weightLimits <- params_belize_analysis$w_a * input$singleSizeLimit ^ params_belize_analysis$w_b
tRes <- input$timeStep
t0 <-tManagement - input$fisheryAge
tVirgin <- t0 - 50
tEnd <- tManagement + input$tProjection


# Virgin Unfished -------------------------------------------------------------------------
# The model of the virgin unfished period
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

## Use regular TidySim here, only applying CC rate to future projections
virginOutput <- tidySim(virginInput)
virginN <- virginOutput$sim@n[dim(virginOutput$sim@n)[1],,]
virginNPP <- virginOutput$sim@n_pp[dim(virginOutput$sim@n_pp)[1],]

# Historic Fishing -------------------------------------------------------------------------
# The model of historic fishing activity
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

historicOutput <- tidySim(historicInput)
baselineN <- historicOutput$sim@n[dim(historicOutput$sim@n)[1],,]
baselineNPP <- historicOutput$sim@n_pp[dim(historicOutput$sim@n_pp)[1],]

# Define Management -------------------------------------------------------------------------
# Create the model parameters of the different management scenarios
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
  
  FInput = list(managementName = "Limit F relative to natural mortality",
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
  
  sizeAndFInput = list(managementName = "Size limits and F limits",
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
# Add seasonal closure -------------------------------------------------------------------------
# Haven't finished the seasonal code yet, add it to the end of the list here

# Run Management Scenarios -------------------------------------------------------------------------
# loops through all management scenarios above and runs the two mizers on them

forwardProjectionInputs2 <- forwardProjectionInputs[as.numeric(1:7)]

## Run non-climate mizer
values$forwardProjectionOutputs <- map(forwardProjectionInputs2,function(x){
  tidySim(x)
  
})

## Run climate mizer
valuesCC$forwardProjectionOutputs <- map(forwardProjectionInputs2,function(x){
    tidySimCC(x)
    
  })


# Reorganize results -------------------------------------------------------------------------
# Arrange results in a way that allows for visualization

## Non-climate mizer
values$combinedSpeciesOutput <- values$forwardProjectionOutputs %>%
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

## Climate mizer
valuesCC$combinedSpeciesOutput <- valuesCC$forwardProjectionOutputs %>%
  map("simResultsSpecies") %>% 
  bind_rows() %>%
  rbind(historicOutput$simResultsSpecies,
        virginOutput$simResultsSpecies) %>%
  filter(!is.na(Year))

valuesCC$combinedAggregatedOutput <- valuesCC$forwardProjectionOutputs %>%
  map("simResultsAggregated") %>% 
  bind_rows() %>%
  rbind(historicOutput$simResultsAggregated,
        virginOutput$simResultsAggregated) %>%
  filter(!is.na(Year))

speciesVirgin <- rownames(virginN)
valuesCC$virginTidy <- virginN %>% as_data_frame()
valuesCC$virginTidy$Species <- speciesVirgin
valuesCC$virginTidy <- valuesCC$virginTidy %>%
  gather(weightClass,virginAbundance,-Species)
valuesCC$virginTidy$weightClass <- as.numeric(valuesCC$virginTidy$weightClass)


simS<-valuesCC$forwardProjectionOutputs %>% map("sim")
timeStep <- input$tProjection +1
simNames <- names(simS)
mapNames <- map(forwardProjectionInputs,function(x){x$managementName}) %>% as.vector()
valuesCC$sizeSpectrum <- map(simNames,function(x){
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
    left_join(valuesCC$virginTidy) %>%
    mutate(relativeAbundance = n/virginAbundance)
}) %>%
  bind_rows()

# Aggregate Plots -------------------------------------------------------------------------
# All 3 types of plots for one projection type on aggregated species

## SWAP OUT VALUES & VALUESCC TO CHOOSE WHICH TO PLOT
plotter <- values
plotter <- valuesCC

t0 <-tManagement - input$fisheryAge
tVirgin <- t0 - 100
tEnd <- tManagement + input$tProjection


bmin <- min(plotter$combinedAggregatedOutput$Biomass)
bmax <- max(plotter$combinedAggregatedOutput$Biomass)
ymin <- min(plotter$combinedAggregatedOutput$Yield)
ymax <- max(plotter$combinedAggregatedOutput$Yield)

## BIOMASS PLOT:
bPlot <- plotter$combinedAggregatedOutput %>%
  mutate(Biomass = (Biomass - bmin) / (bmax - bmin)) %>% ## Convert from grams to kg
  filter(Year >= tVirgin) %>%
  ggplot(aes(x=Year, y=Biomass,color=Management,group=Management)) +
  scale_colour_manual(values =c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit F relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and F limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit F relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and F limits",
                               "No-take zone",
                               "No change", 
                               "Seasonal closure", 
                               "F based on vulnerability basket")) +
  annotate("rect", xmin = -Inf, xmax = t0-1, ymin = -Inf, ymax = Inf, fill = "lightskyblue", alpha = 0.2) +
  annotate("rect", xmin = t0-1, xmax = tManagement - 1, ymin = -Inf, ymax = Inf, fill = "deepskyblue", alpha = 0.2) +
  geom_line(size=1.5) +
  geom_vline(xintercept = tManagement-1) +
  xlab("Year") +
  ylab("Relative biomass") +
  labs(color = "Management intervention") +
  ggtitle("a. Projected relative biomass")  +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        plot.title=element_text(size=16, face="bold"),
        legend.title = element_text(size=12), legend.position = "none")

## YIELD PLOT:
yPlot <- plotter$combinedAggregatedOutput %>%
  filter(Year >= tVirgin) %>%
  mutate(Yield = (Yield - ymin) / (ymax - ymin)) %>% 
  ggplot(aes(x=Year, y=Yield,color=Management,group=Management)) +
  scale_colour_manual(values = c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit F relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and F limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit F relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and F limits",
                               "No-take zone",
                               "No change", 
                               "Seasonal closure", 
                               "F based on vulnerability basket")) +
  annotate("rect", xmin = -Inf, xmax = t0-1, ymin = -Inf, ymax = Inf, fill = "lightskyblue", alpha = 0.2) +
  annotate("rect", xmin = t0-1, xmax = tManagement - 1, ymin = -Inf, ymax = Inf, fill = "deepskyblue", alpha = 0.2) +
  geom_line(size=1.5) +
  geom_vline(xintercept = tManagement-1) +
  xlab("Year") +
  ylab("Relative catch") +
  labs(color = "Management intervention") +
  ggtitle("b. Projected relative catch") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        plot.title=element_text(size=16, face="bold"),
        legend.title = element_text(size=12), legend.position = "none")


## TRADEOFF PLOT:
tradeoffPlot <- plotter$combinedAggregatedOutput %>%
  filter(Management != "Historic fishing" & Management != "Virgin unfished") %>%
  filter(Year == tEnd) %>%
  ggplot(aes(x=((Biomass - min(plotter$combinedAggregatedOutput$Biomass)) / (max(plotter$combinedAggregatedOutput$Biomass) - min(plotter$combinedAggregatedOutput$Biomass))),
             y=((Yield - min(plotter$combinedAggregatedOutput$Yield)) / (max(plotter$combinedAggregatedOutput$Yield) - min(plotter$combinedAggregatedOutput$Yield))),
             color=Management)) +
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit F relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and F limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit F relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and F limits",
                               "No-take zone",
                               "No change", 
                               "Seasonal closure", 
                               "F based on vulnerability basket")) +
  geom_point(size = 5) +
  xlab("Relative biomass") +
  ylab("Relative catch") +
  labs(color = "Management intervention") +
  ggtitle(paste("c. Tradeoff plot of projected relative \nbiomass and catch in",tEnd,sep=" ")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        plot.title=element_text(size=16, face="bold"),
        legend.title = element_text(size=12, face="bold"), legend.text = element_text(size=12))

Agg<-grid.arrange(bPlot,
                  yPlot,
                  tradeoffPlot,
                  ncol=1)
## Only use the below if you want to save the image
ggsave("figs_54_F1_RS12_TL/AggregatePlot2.jpeg", Agg, height=8,width=8, units="in")
print(Agg)
dev.off()

# Species Plots -------------------------------------------------------------------------
# All 3 types of plots for one projection of a single define species
t0 <-tManagement - input$fisheryAge
tVirgin <- t0 - 100
tEnd <- tManagement + input$tProjection

## SWAP OUT VALUES & VALUESCC TO CHOOSE WHICH TO PLOT
plotter <- values
plotter <- valuesCC

## NEED TO FILL IN A SPECIES THAT YOU WANT TO PLOT. Could look at:
## unique(plotter$combinedSpeciesOutput$Species) for the species list

## Loop that plots and saves ALL spp individually
for (i in 1:nsS){
  input$speciesSelect <- unique(plotter$combinedSpeciesOutput$Species)[i]
  file_name = paste("figs_54_F1_RS12_TL/", input$speciesSelect, ".jpeg", sep="")
  
  ## BIOMASS PLOT:
  bPlot<-plotter$combinedSpeciesOutput %>%
    filter(Species == input$species ) %>%
    mutate(Biomass = (Biomass - min(Biomass)) / (max(Biomass) - min(Biomass))) %>% ## Convert to normalized values
    ggplot(aes(x=Year, y=Biomass,color=Management,group=Management)) +
    scale_colour_manual(values=c("Virgin unfished"="black", 
                                 "Historic fishing"="blue",
                                 "No change"="red",
                                 "Limit F relative to natural mortality" = "#009E73",
                                 "Minimum size limit (for each species)"= "slateblue",
                                 "Size limits and F limits"= "yellow",
                                 "No-take zone" = "orange",
                                 "Single F across species" = "purple",
                                 "Single minimum size limit across species" = "pink",
                                 "F based on vulnerability basket" = "green",
                                 "Seasonal closure" = "teal"), 
                        breaks=c("Virgin unfished", "Historic fishing",
                                 "Limit F relative to natural mortality",
                                 "Single F across species",
                                 "Single minimum size limit across species",
                                 "Minimum size limit (for each species)",
                                 "Size limits and F limits",
                                 "No-take zone",
                                 "No change", 
                                 "Seasonal closure", 
                                 "F based on vulnerability basket")) +
    geom_line(size=1.3) +
    geom_vline(xintercept = tManagement-1) +
    ylab("Relative biomass") +
    labs(color = "Management intervention") +
    annotate("rect", xmin = -Inf, xmax = t0-1, ymin = 10 ^ -Inf, ymax = 10 ^ Inf, fill = "lightskyblue", alpha = 0.2) +
    annotate("rect", xmin = t0-1, xmax = tManagement - 1, ymin = 10 ^ -Inf, ymax = 10 ^ Inf, fill = "deepskyblue", alpha = 0.2) +
    ggtitle(paste("a. Projected relative biomass", input$speciesSelect)) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"), 
          plot.title=element_text(size=16, face="bold"),
          legend.title = element_text(size=12), legend.position = "none")
  
  ## YIELD PLOT:
  yPlot<-plotter$combinedSpeciesOutput %>%
    filter(Species == input$speciesSelect) %>%
    mutate(Yield = (Yield - min(Yield)) / (max(Yield) - min(Yield))) %>% ## Convert from grams to kg
    ggplot(aes(x=Year, y=Yield,color=Management,group=Management)) +
    scale_colour_manual(values=c("Virgin unfished"="black", 
                                 "Historic fishing"="blue",
                                 "No change"="red",
                                 "Limit F relative to natural mortality" = "#009E73",
                                 "Minimum size limit (for each species)"= "slateblue",
                                 "Size limits and F limits"= "yellow",
                                 "No-take zone" = "orange",
                                 "Single F across species" = "purple",
                                 "Single minimum size limit across species" = "pink",
                                 "F based on vulnerability basket" = "green",
                                 "Seasonal closure" = "teal"), 
                        breaks=c("Virgin unfished", "Historic fishing",
                                 "Limit F relative to natural mortality",
                                 "Single F across species",
                                 "Single minimum size limit across species",
                                 "Minimum size limit (for each species)",
                                 "Size limits and F limits",
                                 "No-take zone",
                                 "No change", 
                                 "Seasonal closure", 
                                 "F based on vulnerability basket")) +
    annotate("rect", xmin = -Inf, xmax = t0-1, ymin = 10 ^ -Inf, ymax = 10 ^ Inf, fill = "lightskyblue", alpha = 0.2) +
    annotate("rect", xmin = t0-1, xmax = tManagement - 1, ymin = 10 ^ -Inf, ymax = 10 ^ Inf, fill = "deepskyblue", alpha = 0.2) +
    #facet_grid(Management~.) +
    geom_line(size=1.3) +
    geom_vline(xintercept = tManagement-1) +
    ylab("Relative catch") +
    labs(color = "Management intervention") +
    ggtitle("b. Projected relative catch") +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"), 
          plot.title=element_text(size=16, face="bold"),
          legend.title = element_text(size=12), legend.position = "none")
  
  ## TRADEOFF PLOT:
  tradeoffPlot <- plotter$combinedSpeciesOutput %>%
    filter(Management != "Historic fishing" & Management != "Virgin unfished") %>%
    filter(Species == input$speciesSelect) %>%
    filter(Year == tEnd) %>%
    ggplot(aes(x=((Biomass - min(Biomass)) / (max(Biomass) - min(Biomass))),
               y=((Yield - min(Yield)) / (max(Yield) - min(Yield))),
               color=Management))+
    scale_colour_manual(values=c("Virgin unfished"="black", 
                                 "Historic fishing"="blue",
                                 "No change"="red",
                                 "Limit F relative to natural mortality" = "#009E73",
                                 "Minimum size limit (for each species)"= "slateblue",
                                 "Size limits and F limits"= "yellow",
                                 "No-take zone" = "orange",
                                 "Single F across species" = "purple",
                                 "Single minimum size limit across species" = "pink",
                                 "F based on vulnerability basket" = "green",
                                 "Seasonal closure" = "teal"), 
                        breaks=c("Virgin unfished", "Historic fishing",
                                 "Limit F relative to natural mortality",
                                 "Single F across species",
                                 "Single minimum size limit across species",
                                 "Minimum size limit (for each species)",
                                 "Size limits and F limits",
                                 "No-take zone",
                                 "No change", 
                                 "Seasonal closure", 
                                 "F based on vulnerability basket")) +
    geom_point(size = 5) +
    xlab("Relative biomass") +
    ylab("Relative catch") +
    labs(color = "Management intervention") +
    ggtitle(paste("c. Tradeoff plot of projected relative \nbiomass and catch in",tEnd,sep=" ")) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"), 
          plot.title=element_text(size=16, face="bold"),
          legend.title = element_text(size=12, face="bold"), legend.text = element_text(size=12))
  
  
  xx<-grid.arrange(bPlot,
                   yPlot,
                   tradeoffPlot,
                   ncol=1)
  ggsave(file=file_name, xx, height=8,width=8, units="in")
  print(xx)
  
}
dev.off()

# Size spectrum plots -------------------------------------------------------------------------
# The size spectrum plots for one projection type
t0 <-tManagement - input$fisheryAge
tVirgin <- t0 - 100
tEnd <- tManagement + input$tProjection

## SWAP OUT VALUES & VALUESCC TO CHOOSE WHICH TO PLOT
plotter <- values
plotter <- valuesCC

SSP<-plotter$sizeSpectrum %>%
  group_by(Management,weightClass) %>%
  summarize(virginAbundance = sum(virginAbundance),
            fishedAbundance = sum(n),
            relativeAbundance = fishedAbundance/virginAbundance) %>%
  ggplot(aes(x=weightClass, y=relativeAbundance,color=Management)) +
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit F relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and F limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit F relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and F limits",
                               "No-take zone",
                               "No change", 
                               "Seasonal closure", 
                               "F based on vulnerability basket")) +
  geom_line(size=1.5) +
  ylab(paste("Relative abundance in ",tEnd," compared to unfished",sep="")) +
  xlab("Size class (g)") +
  labs(color = "Management intervention") +
  geom_hline(yintercept=1,linetype=2) +
  scale_x_log10() +
  ggtitle("Size spectrum") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        plot.title=element_text(size=16, face="bold"),
        legend.title = element_text(size=12, face="bold"), legend.text = element_text(size=12),
        legend.position = "right")

ggsave("figs_54_F1_RS12_TL/SizeSpectrumPlot.jpeg", SSP, height=6,width=10, units="in")
print(SSP)
dev.off()

# Prep for CC Compare -------------------------------------------------------------------------
# Merge cc and non-cc results in a way they can be plotted simultaneously

compare <- list()
compare$params_belize<-params_belize_analysis

compare$SpeciesCC <- valuesCC$forwardProjectionOutputs %>%
  map("simResultsSpecies") %>% 
  bind_rows() %>%
  filter(!is.na(Year))
compare$Species <- values$forwardProjectionOutputs %>%
  map("simResultsSpecies") %>% 
  bind_rows() %>%
  filter(!is.na(Year))
compare$Species$CC <- "None"
compare$SpeciesCC$CC <- "Climate"
compare$combinedSpeciesOutput <- compare$Species %>%
  rbind(compare$SpeciesCC)

compare$Aggregated <- values$forwardProjectionOutputs %>%
  map("simResultsAggregated") %>% 
  bind_rows() %>%
  filter(!is.na(Year))
compare$AggregatedCC <- valuesCC$forwardProjectionOutputs %>%
  map("simResultsAggregated") %>% 
  bind_rows() %>%
  filter(!is.na(Year))
compare$Aggregated$CC <- "None"
compare$AggregatedCC$CC <- "Climate"
compare$combinedAggregatedOutput <- compare$Aggregated %>%
  rbind(compare$AggregatedCC)

# Compare Aggregate -------------------------------------------------------------------------
# Comparing the aggregate outputs between climate and non-climate models

## Biomass
bPlot<-compare$combinedAggregatedOutput %>%
  mutate(Biomass = (Biomass - min(Biomass)) / (max(Biomass) - min(Biomass))) %>% ## Convert to normalized values
  ggplot(aes(x=Year, y=Biomass,color=Management, lty=CC)) +
  scale_linetype_manual(values = c(None = "dashed", Climate = "solid"), guide = FALSE) +
  scale_colour_manual(values=c("No change"="red",
                               "Limit F relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and F limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink")) +
  geom_line(size=1.3) +
  ylab("Relative biomass") +
  labs(color = "Management intervention") +
  ggtitle("a. Projected relative biomass under extreme CC") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        plot.title=element_text(size=16, face="bold"),
        legend.title = element_blank(), legend.position = "none")

## Yield
yPlot<-compare$combinedAggregatedOutput %>%
  mutate(Yield = (Yield - min(Yield)) / (max(Yield) - min(Yield))) %>% ## Convert to normalized values
  ggplot(aes(x=Year, y=Yield,color=Management,lty=CC)) +
  scale_linetype_manual(values = c(None = "dashed", Climate = "solid"), guide = FALSE) +
  scale_colour_manual(values=c("No change"="red",
                               "Limit F relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and F limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink")) +
  geom_line(size=1.3) +
  ylab("Relative yield") +
  labs(color = "Management intervention") +
  ggtitle("a. Projected relative yield under extreme CC") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        plot.title=element_text(size=16, face="bold"),
        legend.title = element_blank(), legend.position = "bottom")

Comp<-grid.arrange(bPlot,
                  yPlot,
                  ncol=1)
ggsave("figures/ccX_comp_agg.jpeg", Comp, height=8,width=8, units="in")
print(Comp)
dev.off()

# Compare Species -------------------------------------------------------------------------
# Comparing the species outputs between climate and non-climate models

## Select individual species
input$speciesSelect <- unique(compare$combinedSpeciesOutput$Species)[7]
  
## Biomass
bPlot<-compare$combinedSpeciesOutput %>%
  filter(Species == input$speciesSelect ) %>%
  mutate(Biomass = (Biomass - min(Biomass)) / (max(Biomass) - min(Biomass))) %>% ## Convert to normalized values
  ggplot(aes(x=Year, y=Biomass,color=Management, lty=CC)) +
  scale_linetype_manual(values = c(None = "dashed", Climate = "solid"), guide = FALSE) +
  scale_colour_manual(values=c("No change"="red",
                               "Limit F relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and F limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink")) +
  geom_line(size=1.3) +
  ylab("Relative biomass") +
  labs(color = "Management intervention") +
  ggtitle(paste("a. Projected relative biomass under extreme CC for ", input$speciesSelect)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        plot.title=element_text(size=16, face="bold"),
        legend.title = element_text(size=12), legend.position = "none")

## Yield
yPlot<-compare$combinedSpeciesOutput %>%
  filter(Species == input$speciesSelect ) %>%
  mutate(Yield = (Yield - min(Yield)) / (max(Yield) - min(Yield))) %>% ## Convert to normalized values
  ggplot(aes(x=Year, y=Yield,color=Management,lty=CC)) +
  scale_linetype_manual(values = c(None = "dashed", Climate = "solid"), guide = FALSE) +
  scale_colour_manual(values=c("No change"="red",
                               "Limit F relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and F limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink")) +
  geom_line(size=1.3) +
  ylab("Relative yield") +
  labs(color = "Management intervention") +
  ggtitle(paste("a. Projected relative yield under extreme CC for ", input$speciesSelect)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"), 
        plot.title=element_text(size=16, face="bold"),
        legend.title = element_text(size=12), legend.position = "bottom")

CompSpp<-grid.arrange(bPlot,
                   yPlot,
                   ncol=1)
file_name = paste("figures/ccX_comp_", input$speciesSelect, ".jpeg", sep="")
ggsave(file_name, CompSpp, height=8,width=8, units="in")
print(CompSpp)
dev.off()
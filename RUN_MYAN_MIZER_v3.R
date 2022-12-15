# global.r
## Myanmar multispecies finfish study
## Load packages


rm(list = ls())

setwd("C:/Users/kkleisner/Documents/ContractWork/EDF/Upside/Helmsley/Myanmar/myanmar-mizer/Mizer_Myanmar-Clean-kmk/")
## Load data

library(mizer)
library(tidyverse)
source("functions.R")
library(rhandsontable)
library(gridExtra)


params_belize_default <- read_csv("data/myanmar_model_inputs_54_TL.csv",col_names=TRUE)
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


####MAKE AGGREGATE SPECIES PLOTS:

t0 <-tManagement - input$fisheryAge
tVirgin <- t0 - 100
tEnd <- tManagement + input$tProjection


bmin <- min(values$combinedAggregatedOutput$Biomass)
bmax <- max(values$combinedAggregatedOutput$Biomass)
ymin <- min(values$combinedAggregatedOutput$Yield)
ymax <- max(values$combinedAggregatedOutput$Yield)

##BIOMASS PLOT:
bPlot <- values$combinedAggregatedOutput %>%
  mutate(Biomass = (Biomass - bmin) / (bmax - bmin)) %>% ## Convert from grams to kg
  filter(Year >= tVirgin) %>%
  ggplot(aes(x=Year, y=Biomass,color=Management,group=Management)) +
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit catch relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and catch limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "No change", "Limit catch relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "No-take zone",
                               "Seasonal closure", "Size limits and catch limits",
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
#theme(axis.text.y=element_blank())

##YIELD PLOT:
yPlot <- values$combinedAggregatedOutput %>%
  filter(Year >= tVirgin) %>%
  mutate(Yield = (Yield - ymin) / (ymax - ymin)) %>% 
  ggplot(aes(x=Year, y=Yield,color=Management,group=Management)) +
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit catch relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and catch limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "No change", "Limit catch relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "No-take zone",
                               "Seasonal closure", "Size limits and catch limits",
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
#theme(axis.text.y=element_blank())


##TRADEOFF PLOT:
tradeoffPlot <- values$combinedAggregatedOutput %>%
  filter(Management != "Historic Fishing" & Management != "Virgin Unfished") %>%
  filter(Year == tEnd) %>%
  ggplot(aes(x=((Biomass - min(values$combinedAggregatedOutput$Biomass)) / (max(values$combinedAggregatedOutput$Biomass) - min(values$combinedAggregatedOutput$Biomass))),
             y=((Yield - min(values$combinedAggregatedOutput$Yield)) / (max(values$combinedAggregatedOutput$Yield) - min(values$combinedAggregatedOutput$Yield))),
             color=Management)) +
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit catch relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and catch limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit catch relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and catch limits",
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
#theme(axis.text.y=element_blank(),
#     axis.text.x=element_blank())

Agg<-grid.arrange(bPlot,
             yPlot,
             tradeoffPlot,
             ncol=1)
ggsave("figs_54_F1_RS12_TL/AggregatePlot2.jpeg", Agg, height=8,width=8, units="in")
print(Agg)
dev.off()

####MAKE SINGLE SPECIES PLOTS:
t0 <-tManagement - input$fisheryAge
tVirgin <- t0 - 100
tEnd <- tManagement + input$tProjection

#####NEED TO FILL IN A SPECIES THAT YOU WANT TO PLOT. Could look at:
##### unique(values$combinedSpeciesOutput$Species) for the species list

for (i in 1:nsS){
input$speciesSelect<-unique(values$combinedSpeciesOutput$Species)[i]
file_name = paste("figs_54_F1_RS12_TL/", input$speciesSelect, ".jpeg", sep="")

##BIOMASS PLOT:
bPlot<-values$combinedSpeciesOutput %>%
  filter(Species == input$species ) %>%
  mutate(Biomass = (Biomass - min(Biomass)) / (max(Biomass) - min(Biomass))) %>% ## Convert to normalized values
  ggplot(aes(x=Year, y=Biomass,color=Management,group=Management)) +
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit catch relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and catch limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit catch relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and catch limits",
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

##YIELD PLOT:
yPlot<-values$combinedSpeciesOutput %>%
    filter(Species == input$speciesSelect) %>%
    mutate(Yield = (Yield - min(Yield)) / (max(Yield) - min(Yield))) %>% ## Convert from grams to kg
    ggplot(aes(x=Year, y=Yield,color=Management,group=Management)) +
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit catch relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and catch limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit catch relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and catch limits",
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

##TRADEOFF PLOT:
tradeoffPlot <- values$combinedSpeciesOutput %>%
  filter(Management != "Historic Fishing" & Management != "Virgin Unfished") %>%
  filter(Species == input$speciesSelect) %>%
  filter(Year == tEnd) %>%
  ggplot(aes(x=((Biomass - min(Biomass)) / (max(Biomass) - min(Biomass))),
             y=((Yield - min(Yield)) / (max(Yield) - min(Yield))),
             color=Management))+
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit catch relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and catch limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit catch relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and catch limits",
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

####MAKE SIZE SPECTRUM PLOTS:
t0 <-tManagement - input$fisheryAge
tVirgin <- t0 - 100
tEnd <- tManagement + input$tProjection
SSP<-values$sizeSpectrum %>%
  group_by(Management,weightClass) %>%
  summarize(virginAbundance = sum(virginAbundance),
            fishedAbundance = sum(n),
            relativeAbundance = fishedAbundance/virginAbundance) %>%
  ggplot(aes(x=weightClass, y=relativeAbundance,color=Management)) +
  scale_colour_manual(values=c("Virgin unfished"="black", 
                               "Historic fishing"="blue",
                               "No change"="red",
                               "Limit catch relative to natural mortality" = "#009E73",
                               "Minimum size limit (for each species)"= "slateblue",
                               "Size limits and catch limits"= "yellow",
                               "No-take zone" = "orange",
                               "Single F across species" = "purple",
                               "Single minimum size limit across species" = "pink",
                               "F based on vulnerability basket" = "green",
                               "Seasonal closure" = "teal"), 
                      breaks=c("Virgin unfished", "Historic fishing",
                               "Limit catch relative to natural mortality",
                               "Single F across species",
                               "Single minimum size limit across species",
                               "Minimum size limit (for each species)",
                               "Size limits and catch limits",
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

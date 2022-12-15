## Unit conversion
## Default units should be grams and cm. Convert everything else.
unitConversion <- function(Unit,WeightA,WeightB,FL_TL_a,FL_TL_b,SL_TL_a,SL_TL_b){
  pmap(list(Unit,WeightA,WeightB,FL_TL_a,FL_TL_b,SL_TL_a,SL_TL_b),function(x,y,z,xx,yy,zz,xxx){
    if(is.na(x)) return()
    if (grepl("cm",x) & grepl("kg",x)) {
      temp <- c(1000,1)
      temp <- y * temp[1] * temp[2] ^ z}
    if (grepl("cm",x) & grepl("g",x)) {
      temp <- c(1,1)
      temp <- y * temp[1] * temp[2] ^ z}
    if (grepl("mm",x) & grepl("g",x)) {
      temp <- c(1,10)
      temp <- y * temp[1] * temp[2] ^ z}
    if(grepl("fl",x)) temp[1] <- temp[1] * yy ^ z
    if(grepl("sl",x)) temp[1] <- temp[1] * xxx ^ z
    return(temp)
  })
} %>%
  unlist()

vulnerability<-function(vulVec,fL,fM,fH,fVH){
  pmap(list(vulVec),function(v){
    if (v >= 2.2) f <- fVH
    if (v >=2 & v < 2.2) f <- fH
    if(v >= 1.8 & v < 2) f <- fH
    if(v < 1.8) f <- fL
    return(f)
  })
} %>%
  unlist()

FvMFunc<-function(susc){
  pmap(list(susc),function(s){
    if (s >= 2.2) FvM <- 2
    if (s >=2 & s < 2.2) FvM <- 1.5
    if(s >= 1.8 & s < 2) FvM <- 1
    if(s < 1.8) FvM <- 0.5
    return(FvM)
  })
} %>%
  unlist()

## Forward mizer simulation
tidySim <- function(setup){
  managementName <- setup$managementName
  params_data <- setup$params_data
  params_data$sel_func =setup$sel_func
  params_data$knife_edge_size <-setup$knife_edge_size
  timeRes <- setup$res
  interactions <- setup$interactions
  nws <- setup$nws
  nps <- setup$nps
  nss <- setup$nss
  Params <- MizerParams(params_data[1:nss,],no_w=nws,no_w_pp=nps,interaction = interactions) ## initialize MizerParams object
  ## Define natural mortality
  #Params@species_params$M <- exp(-.2107 -.0824*log(Params@species_params$w_inf) + .6757*log(Params@species_params$k_vb+.4627*log(7.2))) %>%
  #setNames(Params@species_params$gear)## approximately natural mortalty, based on Pauly 1980, 10.1093/icesjms/39.2.175
  time <- length(setup$time) -1
  timeDF <- data_frame(Time = seq(0,time),Year =setup$time)
  e0 <- setup$e0
  eDelta <- setup$eDelta
  eDeltaVec <- c(1,(1+eDelta)^(1:time))
  eMatrix <- eDeltaVec %o% e0
  dimnames(eMatrix) <- list(seq(1:length(eDeltaVec)),Params@species_params$gear)
  sim <- mizer::project(Params, effort = eMatrix, t_save = 1, dt = timeRes, t_max = time, initial_n = setup$initN, intial_n_pp = setup$initNPP)
  
  ## extract biomass from mizer simulation object, make tidy
  biomass <- getBiomass(sim) %>%
    as_data_frame() %>%
    mutate(Time = seq(0,nrow(.)-1)) %>%
    gather(Species,Biomass,everything(),-Time)
  ## extract number of individuals from mizer simulation object, make tidy
  n <- getN(sim) %>%
    as_data_frame() %>%
    mutate(Time = seq(0,nrow(.)-1)) %>%
    gather(Species,n,everything(),-Time)
  ## extract yield from mizer simulation object, make tidy
  yield <- getYield(sim) %>%
    as_data_frame() %>%
    mutate(Time = seq(0,nrow(.)-1)) %>%
    gather(Species,Yield,everything(),-Time)
  ## extract effort from mizer simulation object, make tidy
  effort <<- sim@effort %>%
    as_data_frame() %>%
    mutate(Time = seq(0,nrow(.)-1)) %>%
    gather(Species,Effort,everything(),-Time)
  
  ## join all projection dataframes together for species-level results
  simResultsSpecies <- biomass %>%
    left_join(n,by=c("Time","Species")) %>%
    left_join(yield,by=c("Time","Species")) %>%
    left_join(effort,by=c("Time","Species")) %>%
    left_join(timeDF,by="Time") %>%
    dplyr::select(-time) %>%
    mutate(Management = managementName)
  
  simResultsAggregated <- simResultsSpecies %>%
    group_by(Time) %>%
    summarize(Biomass = sum(Biomass),
              Yield = sum(Yield)) %>%
    left_join(timeDF,by="Time") %>%
    dplyr::select(-Time) %>%
    mutate(Management = managementName)
  
  return(list(sim=sim,
              simResultsSpecies=simResultsSpecies,
              simResultsAggregated=simResultsAggregated))
}

#Mizer simulation with climate impact
tidySimCC <- function(setup){
  managementName <- setup$managementName
  params_data <- setup$params_data
  params_data$sel_func =setup$sel_func
  params_data$knife_edge_size <-setup$knife_edge_size
  timeRes <- setup$res
  interactions <- setup$interactions
  nws <- setup$nws
  nps <- setup$nps
  nss <- setup$nss
  Params <- MizerParams(params_data[1:nss,],no_w=nws,no_w_pp=nps,interaction = interactions) ## initialize MizerParams object
  ## Define natural mortality
  #Params@species_params$M <- exp(-.2107 -.0824*log(Params@species_params$w_inf) + .6757*log(Params@species_params$k_vb+.4627*log(7.2))) %>%
  #setNames(Params@species_params$gear)## approximately natural mortalty, based on Pauly 1980, 10.1093/icesjms/39.2.175
  time <- length(setup$time) -1
  timeDF <- data_frame(Time = seq(0,time),Year =setup$time)
  e0 <- setup$e0
  eDelta <- setup$eDelta
  eDeltaVec <- c(1,(1+eDelta)^(1:time))
  eMatrix <- eDeltaVec %o% e0
  dimnames(eMatrix) <- list(seq(1:length(eDeltaVec)),Params@species_params$gear)
  sim <- projectcc(Params, effort = eMatrix, impact = impact, t_save = 1, dt = timeRes, t_max = time, initial_n = setup$initN, intial_n_pp = setup$initNPP)

  ## extract biomass from mizer simulation object, make tidy
  biomass <- getBiomass(sim) %>%
    as_data_frame() %>%
    mutate(Time = seq(0,nrow(.)-1)) %>%
    gather(Species,Biomass,everything(),-Time)
  ## extract number of individuals from mizer simulation object, make tidy
  n <- getN(sim) %>%
    as_data_frame() %>%
    mutate(Time = seq(0,nrow(.)-1)) %>%
    gather(Species,n,everything(),-Time)
  ## extract yield from mizer simulation object, make tidy
  yield <- getYield(sim) %>%
    as_data_frame() %>%
    mutate(Time = seq(0,nrow(.)-1)) %>%
    gather(Species,Yield,everything(),-Time)
  ## extract effort from mizer simulation object, make tidy
  effort <<- sim@effort %>%
    as_data_frame() %>%
    mutate(Time = seq(0,nrow(.)-1)) %>%
    gather(Species,Effort,everything(),-Time)

  ## join all projection dataframes together for species-level results
  simResultsSpecies <- biomass %>%
    left_join(n,by=c("Time","Species")) %>%
    left_join(yield,by=c("Time","Species")) %>%
    left_join(effort,by=c("Time","Species")) %>%
    left_join(timeDF,by="Time") %>%
    dplyr::select(-time) %>%
    mutate(Management = managementName)

  simResultsAggregated <- simResultsSpecies %>%
    group_by(Time) %>%
    summarize(Biomass = sum(Biomass),
              Yield = sum(Yield)) %>%
    left_join(timeDF,by="Time") %>%
    dplyr::select(-Time) %>%
    mutate(Management = managementName)

  return(list(sim=sim,
              simResultsSpecies=simResultsSpecies,
              simResultsAggregated=simResultsAggregated))
}

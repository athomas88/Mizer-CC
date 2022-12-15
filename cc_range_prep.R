library(mizer)
library(tibble)

#Singular rate
df <- read.csv("C:/Users/athomas/Documents/Myanmar/Mizer/data/Myanmar_rangearea_cc.csv")
inputs <- read.csv("C:/Users/athomas/Documents/Myanmar/Mizer/data/myanmar_model_inputs_gm.csv")

yrtot <- 88
df[["cRate"]] <- (df[["yr00"]] / df[["yr12"]]) ^ (1 / yrtot)
df <- df[,c(1,3,24)]
df[["Species"]] <- gsub("_", " ", df[["Species"]])
impact <- df[df$RCP==8.5,]
row.names(impact) <- impact[[1]]


#Year to year variation based on 5 year averages
df <- read.csv("C:/Users/athomas/Documents/Myanmar/Mizer/data/Myanmar_change_cc.csv")
rcp85 <- df[df$RCP==8.5,]
row.names(rcp85) <- rcp85[[1]]
rcp85 <- rcp85[,-c(1,2,3,23,24)]
yr1 <- c(2012,2015,2020,2025,2030,2035,2040,2045,2050,2055,2060,2065,2070,2075,2080,2085,2090,2095,2100)
names(rcp85) <- yr1

nums <- read.csv("C:/Users/athomas/Documents/Myanmar/Mizer/data/Myanmar_rangearea_cc.csv")

for (i in yr1) {
  n  <- i
  i <- as.character(i)
  if (n == 2012){
    rcp85[[i]] <- 0
  }
  else if (n==2015){
    y <- seq(n-2,n-1)
    c <- rcp85[[i]] ^ (1 / (length(y+1)))
    for (e in y){
      e <- as.character(e)
      rcp85[[e]] <- 1 - c
    }
    rcp85[[i]] <- 1 - c
  }
  else {
    y <- seq(n-4,n-1)
    c <- rcp85[[i]] ^ (1 / (length(y+1)))
    for (e in y){
      e <- as.character(e)
      rcp85[[e]] <- 1 - c
    }
    rcp85[[i]] <- 1 - c
  }
}


yrs <- 2012:2100
order <- c()
for (i in yrs){
  order <- c(order, as.character(i))
}
rcp85 <- rcp85[order]


#Mizer code

getCImpact <- function(object, impact, time_range) {
  if (is(object, "MizerSim")) {
    sim <- object
    if (missing(time_range)) {
      time_range <- dimnames(sim@effort)$time
    }
    time_elements <- get_time_elements(sim, time_range, slot_name = "effort")
    f_mort_gear <- getCImpact(sim@params, impact)
    return(f_mort_gear[time_elements, , , , drop = FALSE])
  } else {
    params <- object
    spp <- params@species_params[["species"]]
    impact <- impact[spp,3]
    if (is(impact, "numeric")) {
      no_gear <- dim(params@catchability)[1]
      # If a single value, just repeat it for all gears
      if (length(impact) == 1) {
        impact <- rep(impact, no_gear)
      }
      if (length(impact) != no_gear) {
        stop("impact must be a single value or a vector as long as the number of gears\n")
      }
      # Streamlined for speed increase - note use of recycling
      out <- array(rep(params@catchability, 100), dim = dim(params@selectivity))
      out[] <- impact * c(params@catchability) * c(params@selectivity)
      return(out)
    } else {
      # assuming impact is a matrix, and object is of MizerParams class
      no_gear <- dim(params@catchability)[1]
      if (dim(impact)[2] != no_gear)
        stop("impact array must have a single value or a vector as long as the number of gears for each time step\n")
      # Make the output array - note that we put time as last dimension
      # and then aperm before returning. This is because of the order of
      # the values when we call the other getCImpact function.
      # Fill it up by calling the other function and passing in each line
      # of the impact matrix
      out <- array(NA, dim = c(dim(params@selectivity), dim(impact)[1]),
                   dimnames = c(dimnames(params@selectivity),
                                list(time = dimnames(impact)[[1]])))
      out[] <- apply(impact, 1, function(x) getCImpact(params, x))
      out <- aperm(out, c(4, 1, 2, 3))
      return(out)
    }
  }
}

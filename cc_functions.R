
## Climate Impact

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

## Functions


setGeneric('getCImpact', function(object, impact, ...)
  standardGeneric('getCImpact'))


setMethod('getCImpact', signature(object='MizerParams'),
          function(object, impact, ...){
            spp <- object@species_params[["sciname"]] #Must be 'species' for RUN code, and 'sciname' for Shiny
            impact <- impact[spp,3]
            no_gear <- dim(object@catchability)[1]
            # If a single value, just repeat it for all gears
            if(length(impact) == 1)
              impact <- rep(impact, no_gear)
            if (length(impact) != no_gear)
              stop("impact must be a single value or a vector as long as the number of gears\n")
            # Streamlined for speed increase - note use of recycling
            out <- array(rep(object@catchability, 100), dim = dim(object@selectivity))
            out[] <- impact * c(object@catchability) * c(object@selectivity)
            out <- colSums(out)
            return(out)
          }
)

setGeneric('getZcc', function(object, n, n_pp, effort, impact, m2)
  standardGeneric('getZcc'))

#' \code{getZ} method with \code{m2} argument.

# Called from project()
setMethod('getZcc', signature(object='MizerParams', n = 'matrix', 
                            n_pp = 'numeric', effort='numeric', m2 = 'matrix'),
          function(object, n, n_pp, effort, impact, m2){
            if (!all(dim(m2) == c(nrow(object@species_params),length(object@w)))){
              stop("m2 argument must have dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),")")
            }
            return(m2 + object@mu_b + getFMort(object, effort = effort) + getCImpact(object,impact = impact))
          }
)

#' \code{getZ} method without \code{m2} argument.

setMethod('getZcc', signature(object='MizerParams', n = 'matrix', 
                            n_pp = 'numeric', effort='numeric', m2 = 'missing'),
          function(object, n, n_pp, effort, impact){
            m2 <- getM2(object, n=n, n_pp=n_pp)
            z <- getZcc(object, n=n, n_pp=n_pp, effort=effort, impact=impact, m2=m2)
            return(z)
          }
)

setGeneric('projectcc', function(object, effort, ...)
  standardGeneric('projectcc'))

setMethod('projectcc', signature(object='MizerParams', effort='array'),
          function(object, effort, impact, t_save=1, dt=0.1, initial_n=object@initial_n, 
                   initial_n_pp=object@initial_n_pp, shiny_progress = NULL, ...){
            validObject(object)
            # Check that number and names of gears in effort array is same as in MizerParams object
            no_gears <- dim(object@catchability)[1]
            if(dim(effort)[2] != no_gears){
              no_gears_error_message <- paste("The number of gears in the effort array (length of the second dimension = ", dim(effort)[2], ") does not equal the number of gears in the MizerParams object (", no_gears, ").", sep="")
              stop(no_gears_error_message)
            }
            gear_names <- dimnames(object@catchability)[[1]]
            if(!all(gear_names %in% dimnames(effort)[[2]])){
              gear_names_error_message <- paste("Gear names in the MizerParams object (", paste(gear_names, collapse=", "), ") do not match those in the effort array.", sep="")
              stop(gear_names_error_message)
            }
            # Sort effort array to match order in MizerParams
            effort <- effort[,gear_names, drop=FALSE]
            
            # Blow up time dimension of effort array
            # i.e. effort might have been passed in using time steps of 1, but actual dt = 0.1, so need to blow up
            if (is.null(dimnames(effort)[[1]])){
              stop("The time dimname of the effort argument must be numeric.")
            }
            if (any(is.na(as.numeric(dimnames(effort)[[1]])))){
              stop("The time dimname of the effort argument must be numeric.")
            }
            time_effort <- as.numeric(dimnames(effort)[[1]])
            if (is.unsorted(time_effort)) {
              stop("The time dimname of the effort argument should be increasing.")
            }
            t_max <- time_effort[length(time_effort)]
            # Blow up effort so that rows are dt spaced
            time_effort_dt <- seq(from = time_effort[1], to = t_max, by = dt)
            effort_dt <- t(array(NA, dim = c(length(time_effort_dt), dim(effort)[2]), dimnames=list(time = time_effort_dt, dimnames(effort)[[2]])))
            for (i in 1:(length(time_effort)-1)){
              effort_dt[,time_effort_dt >= time_effort[i]] <- effort[i,]
            }
            effort_dt <- t(effort_dt)
            
            # Make the MizerSim object with the right size
            # We only save every t_save steps
            # Divisibility test needs to be careful about machine rounding errors,
            # see https://github.com/sizespectrum/mizer/pull/2
            if((t_save < dt) || !isTRUE(all.equal((t_save - round(t_save / dt) * dt), 0)))
              stop("t_save must be a positive multiple of dt")
            t_skip <- round(t_save/dt)
            t_dimnames_index <- seq(1, to = length(time_effort_dt), by = t_skip)
            t_dimnames <- time_effort_dt[t_dimnames_index]
            sim <- MizerSim(object, t_dimnames = t_dimnames) 
            # Fill up the effort array
            sim@effort[] <- effort_dt[t_dimnames_index,]
            
            # Set initial population
            sim@n[1,,] <- initial_n 
            sim@n_pp[1,] <- initial_n_pp
            
            # Handy things
            no_sp <- nrow(sim@params@species_params) # number of species
            no_w <- length(sim@params@w) # number of fish size bins
            idx <- 2:no_w
            # If no w_min_idx column in species_params, issue error
            if (!("w_min_idx" %in% names(sim@params@species_params))) {
              stop("w_min_idx column missing in species params")
            }
            # Hacky shortcut to access the correct element of a 2D array using 1D notation
            # This references the egg size bracket for all species, so for example
            # n[w_minidx_array_ref] = n[,w_min_idx]
            w_min_idx_array_ref <- (sim@params@species_params$w_min_idx-1) * no_sp + (1:no_sp)
            
            # sex ratio - DO SOMETHING LATER WITH THIS
            sex_ratio <- 0.5
            
            # Matrices for solver
            A <- matrix(0,nrow=no_sp,ncol=no_w)
            B <- matrix(0,nrow=no_sp,ncol=no_w)
            S <- matrix(0,nrow=no_sp,ncol=no_w)
            
            # initialise n and nPP
            # We want the first time step only but cannot use drop as there may only be a single species
            n <- array(sim@n[1,,],dim=dim(sim@n)[2:3])
            dimnames(n) <- dimnames(sim@n)[2:3]
            n_pp <- sim@n_pp[1,]
            t_steps <- dim(effort_dt)[1]-1
            # Set up progress bar
            pb <- progress::progress_bar$new(
              format = "[:bar] :percent ETA: :eta",
              total = length(t_dimnames_index), width= 60)
            if (hasArg(shiny_progress)) {
              # We have been passed a shiny progress object
              shiny_progress$set(message = "Running simulation", value = 0)
              proginc <- 1/length(t_dimnames_index)
            }
            for (i_time in 1:t_steps){
              # Do it piece by piece to save repeatedly calling methods
              # Calculate amount E_{a,i}(w) of available food
              phi_prey <- getPhiPrey(sim@params, n=n, n_pp=n_pp)
              # Calculate amount f_i(w) of food consumed
              feeding_level <- getFeedingLevel(sim@params, n=n, n_pp=n_pp, phi_prey=phi_prey)
              # Calculate the predation rate
              pred_rate <- getPredRate(sim@params, n=n, n_pp=n_pp, feeding_level=feeding_level)
              # Calculate predation mortality on fish \mu_{p,i}(w)
              #m2 <- getM2(sim@params, n=n, n_pp=n_pp, pred_rate=pred_rate)
              m2 <- getM2(sim@params, pred_rate=pred_rate)
              # Calculate total mortality \mu_i(w)
              z <- getZcc(sim@params, n=n, n_pp=n_pp, effort=effort_dt[i_time,], impact = impact, m2=m2)
              # Calculate predation mortality on the background spectrum
              m2_background <- getM2Background(sim@params, n=n, n_pp=n_pp, pred_rate=pred_rate)
              # Calculate the resources available for reproduction and growth
              e <- getEReproAndGrowth(sim@params, n=n, n_pp=n_pp, feeding_level=feeding_level)
              # Calculate the resources for reproduction
              e_spawning <- getESpawning(sim@params, n=n, n_pp=n_pp, e=e)
              # Calculate the growth rate g_i(w)
              e_growth <- getEGrowth(sim@params, n=n, n_pp=n_pp, e_spawning=e_spawning, e=e)
              # R_{p,i}
              rdi <- getRDI(sim@params, n=n, n_pp=n_pp, e_spawning=e_spawning, sex_ratio=sex_ratio)
              # R_i
              rdd <- getRDD(sim@params, n=n, n_pp=n_pp, rdi=rdi, sex_ratio=sex_ratio)
              
              # Iterate species one time step forward:
              # See Ken's PDF
              # A_{ij} = - g_i(w_{j-1}) / dw_j dt
              A[,idx] <- sweep(-e_growth[,idx-1,drop=FALSE]*dt, 2, sim@params@dw[idx], "/")
              # B_{ij} = 1 + g_i(w_j) / dw_j dt + \mu_i(w_j) dt
              B[,idx] <- 1 + sweep(e_growth[,idx,drop=FALSE]*dt,2,sim@params@dw[idx],"/") + z[,idx,drop=FALSE]*dt
              # S_{ij} <- N_i(w_j)
              S[,idx] <- n[,idx,drop=FALSE]
              # Boundary condition upstream end (recruitment)
              B[w_min_idx_array_ref] <- 1+e_growth[w_min_idx_array_ref]*dt/sim@params@dw[sim@params@species_params$w_min_idx]+z[w_min_idx_array_ref]*dt
              # Update first size group of n
              n[w_min_idx_array_ref] <- (n[w_min_idx_array_ref] + rdd*dt/sim@params@dw[sim@params@species_params$w_min_idx]) / B[w_min_idx_array_ref]
              # Update n
              # for (i in 1:no_sp) # number of species assumed small, so no need to vectorize this loop over species
              #     for (j in (sim@params@species_params$w_min_idx[i]+1):no_w)
              #         n[i,j] <- (S[i,j] - A[i,j]*n[i,j-1]) / B[i,j]
              
              n <- inner_project_loop(no_sp=no_sp, no_w=no_w, n=n, A=A, B=B, S=S,
                                      w_min_idx=sim@params@species_params$w_min_idx)
              
              # Dynamics of background spectrum uses a semi-chemostat model (de Roos - ask Ken)
              # We use the exact solution under the assumption of constant mortality during timestep
              tmp <- (sim@params@rr_pp * sim@params@cc_pp / (sim@params@rr_pp + m2_background))
              n_pp <- tmp - (tmp - n_pp) * exp(-(sim@params@rr_pp+m2_background)*dt)
              
              # Store results only every t_step steps.
              store <- t_dimnames_index %in% (i_time+1)
              if (any(store)){
                # Advance progress bar
                pb$tick()
                if (hasArg(shiny_progress)) {
                  shiny_progress$inc(amount = proginc)
                }
                # Store result
                sim@n[which(store),,] <- n 
                sim@n_pp[which(store),] <- n_pp
              }
            }
            # and end
            return(sim)
          }
)


inner_project_loop <- function(no_sp, no_w, n, A, B, S, w_min_idx) {
  .Call('_mizer_inner_project_loop', PACKAGE = 'mizer', no_sp, no_w, n, A, B, S, w_min_idx)
}
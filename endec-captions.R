#in biting_process.R
endec_adjusted_mortality<-function(parameters, species, timestep){
  if (parameters$endec){
    # if (timestep == 130){browser()}
    if (timestep %in% parameters$endec_on){
      eff_len <- timestep-max(parameters$endec_ts[parameters$endec_ts <= timestep]) 
      #eff_len <- seq(1, 23, 1)
      endec_killing <- parameters$mum[[species]] + (parameters$mu_endec[[species]]*exp(-parameters$wane_endec[[species]]*eff_len)) #t - endec_on_t
      #endec_killing <- parameters$mum + (parameters$mu_endec*exp(-parameters$wane_endec*eff_len)) #t - endec_on_t
      return(endec_killing)
    }
    #endec_mu <- ifelse(parameters$endec_coverages[matches] > 0, endec_killing, 0)
    #endec_killing <- parameters$mu_endec[[species]] this works as a switch
    #if (any(matches)){ #add and coverage greater than 0 here: so if the coverage is 0, drop back to the baseline mu
    
    #mu_endec<-(parameters$mum[[species]] + parameters$mu_endec[[species]])*parameters$Q0[[species]]*parameters$endec_coverages[matches] #slightly different to atsb
    
    #mu_atsb<-parameters$mu_atsb[[species]]*parameters$atsb_coverages[matches]
    #return(mu+mu_atsb)
    
    #}
    else {return(parameters$mum[[species]])}
  } 
}

#in mosquito_biology.R
death_rate <- function(f, W, Z, species, parameters,timestep) {
  mum <- parameters$mum[[species]]
  mum <- endec_adjusted_mortality(mum, parameters, species, timestep)
  p1_0 <- exp(-mum * parameters$foraging_time[[species]])
  gonotrophic_cycle <- get_gonotrophic_cycle(species, parameters)
  p2 <- exp(-mum * gonotrophic_cycle) #change what goes into background mort
  p1 <- p1_0 * W / (1 - Z * p1_0) #nets come in through W and Z to affect the mortality rate
  -f * log(p1 * p2)
}

#in vector_control_parameters.R
set_endectocide <- function(
    parameters,
    timesteps,
    coverages, 
    endec_on
) {
  if (length(timesteps) != length(coverages)) {
    stop('timesteps and coverages must align')
  }
  endec_min_age <- 5
  endec_max_age <- 90
  parameters$endec <- TRUE
  parameters$endec_timesteps <- timesteps
  parameters$endec_on <- endec_on
  #parameters$endec_coverages <- coverages*(exp(-endec_min_age/21) - exp(-endec_max_age/21)) #may need to change this with Anna's work
  parameters$endec_coverages <- coverages #just operating it as a switch, background work happens in odin
  parameters
}

#in parameters. R
get_parameters <- function(overrides = list()) {
  parameters <- list(
    dd    = 5,
    dt    = 5,
    da    = 195,
    du    = 110,
    del   = 6.64,
    dl    = 3.72,
    dpl   = .643,
    mup   = .249,
    mum   = .132,
    mu_endec = c(0.2, 0.2, 0.2), #this is the mu_h. Need to check if varies with ivm_cov and Q0.
    wane_endec = c(0.04, 0.04, 0.04),
    sigma_squared   = 1.67,
    n_heterogeneity_groups = 5,
    # immunity decay rates
   # ...

   )}
#' @title Parameterise a bed net strategy
#'
#' @description The model will distribute bed nets at `timesteps` to a random
#' sample of the entire human population. The sample size will be a proportion
#' of the human population taken from the corresponding `coverages`.
#' The sample _can_ contain humans who already have bed nets.
#'
#' All of the sample "use" their bed nets on the timestep after they are
#' distributed. Incomplete usage is not part of this model.
#'
#' If a human in the sample already has a bed net, their bed net will be replaced
#' by a new one.
#'
#' Bed nets will be randomly removed each timestep with a rate of `1 -
#' exp(-1/retention)`
#'
#' The structure for the bed net model is documented in the
#' S.I. of 10.1038/s41467-018-07357-w
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to distribute bed nets
#' @param coverages the proportion of the human population who receive bed nets
#' @param retention the average number of timesteps a net is kept for
#' @param dn0 a matrix of death probabilities for each species over time.
#' With nrows=length(timesteps), ncols=length(species)
#' @param rn a matrix of repelling probabilities for each species over time
#' With nrows=length(timesteps), ncols=length(species)
#' @param rnm a matrix of minimum repelling probabilities for each species over time
#' With nrows=length(timesteps), ncols=length(species)
#' @param gamman a vector of bednet half-lives for each distribution timestep
#' @export
set_bednets <- function(
    parameters,
    timesteps,
    coverages,
    retention,
    dn0,
    rn,
    rnm,
    gamman
) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  lengths <- vnapply(list(coverages, gamman), length)
  if (!all(lengths == length(timesteps))) {
    stop('timesteps and time-varying parameters must align')
  }
  for (x in list(dn0, rn, rnm)) {
    if (ncol(x) != length(parameters$species)) {
      stop('death and repelling probabilities rows need to align with species')
    }
    if (nrow(x) != length(timesteps)) {
      stop('death and repelling probabilities columns need to align with timesteps')
    }
  }
  parameters$bednets <- TRUE
  parameters$bednet_timesteps <- timesteps
  parameters$bednet_coverages <- coverages
  parameters$bednet_dn0 <- dn0
  parameters$bednet_rn <- rn
  parameters$bednet_rnm <- rnm
  parameters$bednet_gamman <- gamman
  parameters$bednet_retention <- retention
  parameters
}

#' @title Parameterise an indoor spraying strategy
#'
#' @description The model will apply indoor spraying at `timesteps` to a random
#' sample of the entire human population. The sample size will be a proportion
#' of the human population taken from the corresponding `coverages`.
#' The sample _can_ contain humans who have already benefited from spraying.
#'
#' If a human in the sample lives in a sprayed house, the efficacy of the
#' spraying will be returned to the maximum.
#'
#' The structure for the indoor residual spraying model is documented in the
#' S.I. of 10.1038/s41467-018-07357-w
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to spray
#' @param coverages the proportion of the population who get indoor
#' spraying
#' @param ls_theta matrix of mortality parameters
#' With nrows=length(timesteps), ncols=length(species)
#' @param ls_gamma matrix of mortality parameters per timestep
#' With nrows=length(timesteps), ncols=length(species)
#' @param ks_theta matrix of feeding success parameters per timestep 
#' With nrows=length(timesteps), ncols=length(species)
#' @param ks_gamma matrix of feeding success parameters per timestep 
#' With nrows=length(timesteps), ncols=length(species)
#' @param ms_theta matrix of deterrence parameters per timestep 
#' With nrows=length(timesteps), ncols=length(species)
#' @param ms_gamma matrix of deterrence parameters per timestep 
#' With nrows=length(timesteps), ncols=length(species)
#' @export
set_spraying <- function(
    parameters,
    timesteps,
    coverages,
    ls_theta,
    ls_gamma,
    ks_theta,
    ks_gamma,
    ms_theta,
    ms_gamma
) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  if (length(coverages) != length(timesteps)) {
    stop('coverages and timesteps must must align')
  }
  decays <- list(
    ls_theta,
    ls_gamma,
    ks_theta,
    ks_gamma,
    ms_theta,
    ms_gamma
  )
  for (x in decays) {
    if (ncol(x) != length(parameters$species)) {
      stop('theta and gamma rows need to align with species')
    }
    if (nrow(x) != length(timesteps)) {
      stop('theta and gamma cols need to align with timesteps')
    }
  }
  parameters$spraying <- TRUE
  parameters$spraying_timesteps <- timesteps
  parameters$spraying_coverages <- coverages
  parameters$spraying_ls_theta <- ls_theta
  parameters$spraying_ls_gamma <- ls_gamma
  parameters$spraying_ks_theta <- ks_theta
  parameters$spraying_ks_gamma <- ks_gamma
  parameters$spraying_ms_theta <- ms_theta
  parameters$spraying_ms_gamma <- ms_gamma
  parameters
}

#' @title Parameterise rescaling of the baseline carrying capacity
#' 
#' @description
#' For a number of reasons we may wish to increase or decrease the carrying capacity
#' for a mosquito species throughout the simulation. For example:
#' 
#' \strong{Larval source management}.
#' 
#' To simulate larval source management, we may wish to decrease the
#' baseline carrying capacity. For example, to simulate a 20% reduction in the 
#' baseline carrying capacity (which may approximate a 20% reduction in larval
#' breeding sites) we can set the scaler to 0.8 at the time of implementation.
#' 
#' \strong{Invasive species}
#' 
#' To capture the expansion of an invasive species into a new niche, 
#' we may wish to increase the baseline carrying capacity for that species at the
#' point of invasion. In this instance we set the scaler to be > 1. Note, that the
#' invading species proportion must be initialised to >0, so in this case set the 
#' initial proportion at a very low value (e.g 0.005), then scale up the carrying
#' capacity at the point of invasion by setting the scaler for that species >>1.
#' 
#' @param parameters a list of parameters to modify 
#' @param timesteps the timesteps when rescaling of baseline carrying capacity will occur
#' @param scaler A matrix of scaling factors for baseline carrying capacity where cols must
#' equal the number of species and rows must align with timesteps
#'
#' @export
set_flexible_carrying_capacity <- function(
    parameters,
    carrying_capacity
){
  if(ncol(carrying_capacity) != length(parameters$species)){
    stop("carrying_capacity cols need to align with number of mosquito species")
  }
  if(nrow(carrying_capacity) != timesteps){
    stop("carrying_capacity rows need to align with timesteps")
  }
  stopifnot(min(carrying_capacity) >= 0)
  stopifnot(min(timesteps) >= 0)
  
  parameters$flexible_carrying_capacity <- TRUE
  parameters$carrying_capacity <- carrying_capacity
  parameters
}
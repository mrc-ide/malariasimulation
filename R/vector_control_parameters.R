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

#' @title Parameterise larval source management
#' 
#' @description Reduces the baseline carrying capacity for a given species by
#' (1 - coverage).
#' 
#' @param parameters the model parameters
#' @param timesteps vector of timesteps for each coverage change
#' @param coverages matrix of coverages 
#' With nrows = length(timesteps), ncols=length(species)
#' 
set_larval_source_management <- function(
    parameters,
    timesteps,
    coverages
){
  stopifnot(nrow(coverages) == length(timesteps))
  stopifnot(ncol(coverages) == length(parameters$species))
  stopifnot(min(timesteps) > 0)
  stopifnot(max(coverages) <= 1)
  parameters$larval_source_management <- TRUE
  parameters$lsm_timesteps <- timesteps
  parameters$lsm_coverages <- coverages
  parameters
}

#' @title Rescale baseline carrying capacity
#' 
#' @description Rescales the baseline carrying capacity for a given species.
#' 
#' @param parameters the model parameters
#' @param timesteps vector of timesteps for each rescale change
#' @param scalers matrix of scalers 
#' With nrows = length(timesteps), ncols=length(species)
#' 
set_rescaled_carrying_capacity <- function(
    parameters,
    timesteps,
    scalers
){
  stopifnot(nrow(scalers) == length(timesteps))
  stopifnot(ncol(scalers) == length(parameters$species))
  stopifnot(min(timesteps) > 0)
  stopifnot(min(scalers) >= 0)
  parameters$rescale_carrying_capacity <- TRUE
  parameters$rcc_timesteps <- timesteps
  parameters$rcc_scalers <- scalers
  parameters
}

#' @title Parameterise custom baseline carrying capacity
#' 
#' @description Allows the user to set a completely flexible and custom
#' carrying capacity for each species
#' 
#' @param parameters the model parameters
#' @param timesteps vector of timesteps for each rescale change
#' @param carrying_capacity matrix of baseline carrying_capacity for each species 
#' With nrows = length(timesteps), ncols=length(species)
#' 
set_flexible_carrying_capacity <- function(
    parameters,
    timesteps,
    carrying_capacity
){
  stopifnot(nrow(carrying_capacity) == length(timesteps))
  stopifnot(ncol(carrying_capacity) == length(parameters$species))
  stopifnot(min(timesteps) > 0)
  stopifnot(min(carrying_capacity) >= 0)

  parameters$flexible_carrying_capacity <- TRUE
  parameters$fcc_timesteps <- timesteps
  parameters$fcc <- carrying_capacity
  parameters
}

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

#' @title Parameterise a habitat management strategy
#'
#' @description The model will simulate larval source habitat management at 
#' `timesteps` to the entire human population. 
#' The impact will reduce adult mosquitoes recruited from pupae to the susceptible
#' cohort.
#'
#' The rapidity of the impact is captured by the lsm_factor function
#' by changing the parameters lsm_rate and deprec_param
#'
#' The amount of impact in terms of reduced adult mosquito densities
#' that the habitat management has is determined by larvi_min which is
#' set to 1 as default when there is no impact.
#'
#' Habitat management effects will be randomly removed each timestep with a rate of `1 -
#' exp(-1/habitat_management_waning)`
#'
#' The structure for the habitat model is documented in the
#' working draft Sherrard-Smith et al (larviciding in Kenya)
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to distribute lsm rounds
#' @param larvi_min the proportion of the mosquito population reduced (0 to 1), 1 is no impact
#' @param lsm_rate the duration til the impact is acheived
#' @param deprec_param the speed on impact taking mosquito densities down to new-normal 
#' @param habitat_management_waning the average number of timesteps LSM works for
#' With nrows=length(timesteps), ncols=length(species)
#' @param larvi_min a matrix of reduction probabilities for each species over time
#' With nrows=length(timesteps), ncols=length(species)
#' @param lsm_rate a matrix of reduction probabilities for each species over time
#' With nrows=length(timesteps), ncols=length(species)
#' @param deprec_param a matrix of reduction probabilities for each species over time
#' @export
set_habitat_management <- function(
  parameters,
  timesteps,
  larvi_min,
  lsm_rate,
  deprec_param
  ) {
  for (x in list(larvi_min, lsm_rate, deprec_param)) {
    if (ncol(x) != length(parameters$species)) {
      stop('habitat management probabilities rows need to align with species')
    }
    if (nrow(x) != length(1)) {
      stop('habitat management probabilities columns needs to be 1')
    }
  }
  parameters$habitat_management <- TRUE
  parameters$habitat_management_timesteps <- timesteps
  parameters$larvi_min <- larvi_min
  parameters$lsm_rate <- lsm_rate
  parameters$deprec_param <- deprec_param
  parameters
}

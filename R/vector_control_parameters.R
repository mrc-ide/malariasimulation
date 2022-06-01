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

#' @title Parameterise a larval source management strategy
#'
#' @description The model will simulate larval source management at 
#' `timesteps` to the entire human population. 
#' The impact will reduce adult mosquitoes recruited from pupae to the susceptible
#' cohort.
#'
#' @param parameters a list of parameters to modify
#' @param habitat_management_timesteps the timesteps at which to distribute lsm rounds
#' @param lsm_new_eqm controls the level of recruitment to the adult mosquito susceptible cohort; 
#' Default = 1, this is no reduction, while 0 simulates mosquito elimination.
#' This must be a matrix of reduction probabilities for each species over time
#' with nrows=length(timesteps), ncols=length(species)
#' @param lsm_rate_alpha a matrix for each species over timesteps
#' With nrows=length(timesteps), ncols=length(species); 
#' Default = -4, highly recommend restricting this between -6 and 5
#' The higher the value the longer it takes for lsm to have an impact
#' @param lsm_rate_beta a matrix for each species over timesteps
#' With nrows=length(timesteps), ncols=length(species); 
#' Default = 0.1, this together with lsm_rate_alpha = -4 will bring about the new equilibrium in about a month, 
#' recommended range 0.5 (prolonged time until new equilibrium level of adult recruitment) 
#' to 0.01 (very rapid decline in mosquito recruitment)
#' New equilibrium defined by lsm_new_eqm parameter estimate
#' @export
set_habitat_management <- function(
  parameters,
  habitat_management_timesteps,
  lsm_new_eqm,
  lsm_rate_alpha,
  lsm_rate_beta
  ) {
  for (x in list(lsm_new_eqm, lsm_rate_alpha, lsm_rate_beta)) {
    if (ncol(x) != length(parameters$species)) {
      stop('habitat management probabilities columns need to align with species')
    }
    if (nrow(x) != length(1)) {
      stop('habitat management probabilities need to have just one row. 
           Only one change in recruitment can be made currently. This corresponds to the time when larval source management is implemented.')
    }
  }
  ## Need to add check for parameters$individual_mosquitoes = FALSE?
  parameters$habitat_management <- TRUE
  parameters$habitat_management_timesteps <- habitat_management_timesteps
  parameters$lsm_new_eqm <- lsm_new_eqm
  parameters$lsm_rate_alpha <- lsm_rate_alpha
  parameters$lsm_rate_beta <- lsm_rate_beta
  parameters
}

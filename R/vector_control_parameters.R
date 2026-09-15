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

#' @title Parameterise an indoor spatial emanator strategy
#'
#' @description The model will apply indoor spatial emanators 
#' at `timesteps` to a random sample of the entire human population. 
#' The sample size will be a proportion of the human population taken 
#' from the corresponding `coverages`.
#' The sample _can_ contain humans who have already benefited from spatial emanators.
#'
#' If a human in the sample lives in a treated house, the efficacy of the
#' spatial emanator will be returned to the maximum impact.
#'
#' The structure for the spatial emanator model will be documented in
#' Kuipou et al (PI Ellie Sherrard-Smith)
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to spray
#' @param coverages the proportion of the population who get emanators
#' @param R0 maximum indoor mosquito repellence at t = 0: matrix with nrows=length(timesteps), ncols=length(species)
#' @param K0 maximum indoor mosquito killing before feeding at t = 0: matrix with nrows=length(timesteps), ncols=length(species)
#' @param M0 maximum indoor mosquito killing after feeding at t = 0: matrix with nrows=length(timesteps), ncols=length(species)
#' @param delta0 maximum outdoor repellence at t = 0: matrix with nrows=length(timesteps), ncols=length(species)
#' @param nR Hill function shape parameter for indoor mosquito repellence: matrix with nrows=length(timesteps), ncols=length(species)
#' @param nK Hill function shape parameter for indoor mosquito killing before feeding: matrix with nrows=length(timesteps), ncols=length(species)
#' @param nM Hill function shape parameter for indoor mosquito killing after feeding: matrix with nrows=length(timesteps), ncols=length(species)
#' @param nD Hill function shape parameter for outdoor mosquito repellence, ncols=length(species)
#' @param eps_R Hill function probability of indoor mosquito repellence at T_ref (equivalent to T_50, when T_ref = 50): matrix with nrows=length(timesteps), ncols=length(species)
#' @param eps_K Hill function probability of indoor mosquito killing before feeding at T_ref (equivalent to T_50, when T_ref = 50): matrix with nrows=length(timesteps), ncols=length(species)
#' @param eps_M Hill function probability of indoor mosquito killing after feeding at T_ref (equivalent to T_50, when T_ref = 50): matrix with nrows=length(timesteps), ncols=length(species)
#' @param eps_d Hill function probability of outdoor mosquito repellence at T_ref (equivalent to T_50, when T_ref = 50): matrix with nrows=length(timesteps), ncols=length(species)
#' @param T_ref day at which eps parameters (R, K, M, det) were estimated: matrix with nrows=length(timesteps), ncols=length(species)
#' @export
set_spatial_emanator <- function(
    parameters,
    timesteps,
    coverages,
    R0, nR, eps_R,
    K0, nK, eps_K,
    M0, nM, eps_M,
    delta0, nd, eps_d,
    T_ref
) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  if (length(coverages) != length(timesteps)) {
    stop('coverages and timesteps must must align')
  }
  
  decay_parameters <- list(
    R0, nR, eps_R,
    K0, nK, eps_K,
    M0, nM, eps_M,
    delta0, nd, eps_d,
    T_ref
  )

  for (x in decay_parameters) {
    if (ncol(x) != length(parameters$species)) {
      stop('theta and gamma rows need to align with species')
    }
    if (nrow(x) != length(timesteps)) {
      stop('theta and gamma cols need to align with timesteps')
    }
  }
  
  if (length(coverages) != length(timesteps)) {
    stop('coverages and timesteps must must align')
  }
  
  parameters$spatial_emanator <- TRUE
  parameters$spatial_emanator_timesteps <- timesteps
  parameters$spatial_emanator_coverages <- coverages
  parameters$spatial_emanator_R0 <- R0
  parameters$spatial_emanator_K0 <- K0
  parameters$spatial_emanator_M0 <- M0
  parameters$spatial_emanator_delta0 <- delta0
  parameters$spatial_emanator_nR <- nR
  parameters$spatial_emanator_nK <- nK
  parameters$spatial_emanator_nM <- nM
  parameters$spatial_emanator_nd <- nd
  parameters$spatial_emanator_eps_R <- eps_R
  parameters$spatial_emanator_eps_K <- eps_K
  parameters$spatial_emanator_eps_M <- eps_M
  parameters$spatial_emanator_eps_d <- eps_d
  parameters$spatial_emanator_T_ref <- T_ref
  
  parameters
}

#' @title Parameterise an outdoor spatial emanator strategy
#'
#' @description The model will apply outdoor spatial emanators 
#' at `timesteps` to a random sample of the entire human population. 
#' The sample size will be a proportion of the human population taken 
#' from the corresponding `coverages`.
#' The sample _can_ contain humans who have already benefited from spatial emanators.
#'
#' If a human in the sample lives in a treated house, the efficacy of the
#' spatial emanator will be returned to the maximum impact.
#'
#' The structure for the spatial emanator model will be documented in
#' Kuipou et al (PI Ellie Sherrard-Smith)
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to spray
#' @param coverages the proportion of the population who get emanators
#' @param dse_out_theta matrix of mortality impact outdoor parameters
#' With nrows=length(timesteps), ncols=length(species)
#' @param dse_out_gamma matrix of mortality impact outdoor parameters per timestep
#' With nrows=length(timesteps), ncols=length(species)
#' @param rse_out_theta matrix of repellence outdoor parameters
#' With nrows=length(timesteps), ncols=length(species)
#' @param rse_out_gamma matrix of repellence outdoor parameters per timestep
#' With nrows=length(timesteps), ncols=length(species)
#' @export
set_spatial_emanator_outdoor <- function(
    parameters,
    timesteps,
    coverages,
    rse_out_theta,
    rse_out_gamma,
    dse_out_theta,
    dse_out_gamma
) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  if (length(coverages) != length(timesteps)) {
    stop('coverages and timesteps must must align')
  }
  decays <- list(
    rse_out_theta,
    rse_out_gamma,
    dse_out_theta,
    dse_out_gamma
  )
  for (x in decays) {
    if (ncol(x) != length(parameters$species)) {
      stop('theta and gamma rows need to align with species')
    }
    if (nrow(x) != length(timesteps)) {
      stop('theta and gamma cols need to align with timesteps')
    }
  }
  parameters$spatial_emanator_outdoor <- TRUE
  parameters$spatial_emanator_outdoor_timesteps <- timesteps
  parameters$spatial_emanator_outdoor_coverages <- coverages
  parameters$spatial_emanator_out_theta <- rse_out_theta
  parameters$spatial_emanator_out_gamma <- rse_out_gamma
  parameters$spatial_emanator_mort_out_theta <- dse_out_theta
  parameters$spatial_emanator_mort_out_gamma <- dse_out_gamma
  
  parameters
}

#' @title Parameterise a semiochemical strategy
#'
#' @description The model will deploy semiochemicals at `timesteps` to a random
#' sample of the entire vector population. 
#'
#' Theory and data are courtesey of Noushin Emami (LSTM in press)
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to deploy semiochemicals
#' @param semiochemical_effect matrix of impact on blood feeding rates, 
#' With nrows=length(timesteps), ncols=length(species)
#' @export
set_semiochemical <- function(
    parameters,
    semiochemical_effect,
    timesteps
) {
  if (nrow(semiochemical_effect) != length(timesteps)) {
    stop('semiochemical_effect and timesteps must align')
  }
  parameters$semiochemical <- TRUE
  parameters$semiochemical_timesteps <- timesteps
  parameters$semiochemical_effect <- semiochemical_effect
  parameters
}

#' @title Parameterise custom baseline carrying capacity
#' 
#' @description Allows the user to set a completely flexible and custom
#' carrying capacity for each species
#' 
#' @param parameters the model parameters
#' @param timesteps vector of timesteps for each rescale change
#' @param carrying_capacity_scalers matrix of scaling factors to scale the baseline 
#' carrying capacity for each species with nrows = length(timesteps),
#'  ncols = length(species)
#' 
#' @export
set_carrying_capacity <- function(
    parameters,
    timesteps,
    carrying_capacity_scalers
){
  stopifnot(nrow(carrying_capacity_scalers) == length(timesteps))
  stopifnot(ncol(carrying_capacity_scalers) == length(parameters$species))
  stopifnot(min(timesteps) > 0)
  stopifnot(min(carrying_capacity_scalers) >= 0)
  
  parameters$carrying_capacity <- TRUE
  parameters$carrying_capacity_timesteps <- timesteps
  parameters$carrying_capacity_scalers <- carrying_capacity_scalers
  parameters
}

#' Get initialised carrying capacity for each species
#'
#' @param parameters the model parameters
#'
#' @return a vector of carrying initialised carrying capacity estimates for
#' each vector species 
#' @export
get_init_carrying_capacity <- function(parameters){
  init_cc <- sapply(1:length(parameters$species), function(x){
    p <- parameters$species_proportions[[x]]
    m <- p * parameters$total_M
    calculate_carrying_capacity(parameters, m, x)
  })
  names(init_cc) <- parameters$species
  return(init_cc)
}

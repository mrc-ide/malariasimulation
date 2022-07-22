#' @title Parameterise a Mass Drug Administration
#' @param parameters a list of parameters to modify
#' @param drug the index of the drug to administer
#' @param timesteps a vector of timesteps for each round of mda
#' @param coverages a vector of the proportion of the target population who receive each
#' round
#' @param min_ages a vector of minimum ages of the target population for each round exclusive (in timesteps)
#' @param max_ages a vector of maximum ages of the target population for each round exclusive (in timesteps)
#' drug
#' @export
set_mda <- function(
  parameters,
  drug,
  timesteps,
  coverages,
  min_ages,
  max_ages
  ) {
  
  if(length(coverages) != length(timesteps)){
    stop("coverages and timesteps do no align")
  }
  if(length(min_ages) != length(timesteps)){
    stop("minimum ages and timesteps do no align")
  }
  if(length(max_ages) != length(timesteps)){
    stop("maximum ages and timesteps do no align")
  }
  
  parameters$mda <- TRUE
  parameters$mda_drug <- drug
  parameters$mda_timesteps <- timesteps
  parameters$mda_coverages <- coverages
  parameters$mda_min_ages <- min_ages
  parameters$mda_max_ages <- max_ages
  parameters
}

#' @title Parameterise a Seasonal Malaria Chemoprevention
#' @param parameters a list of parameters to modify
#' @param drug the index of the drug to administer
#' @param timesteps a vector of timesteps for each round of smc
#' @param coverages a vector of the proportion of the target population who receive each
#' round
#' @param min_ages a vector of minimum ages of the target population for each round exclusive (in timesteps)
#' @param max_ages a vector of maximum ages of the target population for each round exclusive (in timesteps)
#' drug
#' @export
set_smc <- function(
  parameters,
  drug,
  timesteps,
  coverages,
  min_ages,
  max_ages
  ) {
  
  if(length(coverages) != length(timesteps)){
    stop("coverages and timesteps do no align")
  }
  if(length(min_ages) != length(timesteps)){
    stop("minimum ages and timesteps do no align")
  }
  if(length(max_ages) != length(timesteps)){
    stop("maximum ages and timesteps do no align")
  }
  
  parameters$smc <- TRUE
  parameters$smc_drug <- drug
  parameters$smc_timesteps <- timesteps
  parameters$smc_coverages <- coverages
  parameters$smc_min_ages <- min_ages
  parameters$smc_max_ages <- max_ages
  parameters
}

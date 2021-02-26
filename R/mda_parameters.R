#' @title Parameterise a Mass Drug Administration
#' @param parameters a list of parameters to modify
#' @param drug the index of the drug to administer
#' @param timesteps a vector of timesteps for each round of mda
#' @param coverages the proportion of the target population who recieve each
#' round
#' @param min_age the minimum age of the target population exclusive (in timesteps)
#' @param max_age the maximum age of the target population exclusive (in timesteps)
#' drug
#' @export
set_mda <- function(
  parameters,
  drug,
  timesteps,
  coverages,
  min_age,
  max_age
  ) {
  parameters$mda <- TRUE
  parameters$mda_drug <- drug
  parameters$mda_timesteps <- timesteps
  parameters$mda_coverages <- coverages
  parameters$mda_min_age <- min_age
  parameters$mda_max_age <- max_age
  parameters
}

#' @title Parameterise a Seasonal Malaria Chemoprevention
#' @param parameters a list of parameters to modify
#' @param drug the index of the drug to administer
#' @param timesteps a vector of timesteps for each round of smc
#' @param coverages the proportion of the target population who recieve each
#' round
#' @param min_age the minimum age of the target population exclusive (in timesteps)
#' @param max_age the maximum age of the target population exclusive (in timesteps)
#' drug
#' @export
set_smc <- function(
  parameters,
  drug,
  timesteps,
  coverages,
  min_age,
  max_age
  ) {
  parameters$smc <- TRUE
  parameters$smc_drug <- drug
  parameters$mda_timesteps <- timesteps
  parameters$smc_coverages <- coverages
  parameters$smc_min_age <- min_age
  parameters$smc_max_age <- max_age
  parameters
}

#' @title Parameterise a Mass Drug Administration
#' @param parameters a list of parameters to modify
#' @param drug the index of the drug to administer
#' @param start the timestep to start
#' @param end the last timestep for the intervention
#' @param frequency the number of timsteps between doses
#' @param min_age the minimum age of the target population exclusive (in timesteps)
#' @param max_age the maximum age of the target population exclusive (in timesteps)
#' @param coverage the proportion of the target population who will recieve the
#' drug
#' @export
set_mda <- function(
  parameters,
  drug,
  start,
  end,
  frequency,
  min_age,
  max_age,
  coverage
  ) {
  parameters$mda <- TRUE
  parameters$mda_drug <- drug
  parameters$mda_start <- start
  parameters$mda_end <- end
  parameters$mda_frequency <- frequency
  parameters$mda_min_age <- min_age
  parameters$mda_max_age <- max_age
  parameters$mda_coverage <- coverage
  parameters
}

#' @title Parameterise a Seasonal Malaria Chemoprevention
#' @param parameters a list of parameters to modify
#' @param drug the index of the drug to administer
#' @param start the timestep to start
#' @param end the last timestep for the intervention
#' @param frequency the number of timsteps between doses
#' @param min_age the minimum age of the target population exclusive (in timesteps)
#' @param max_age the maximum age of the target population exclusive (in timesteps)
#' @param coverage the proportion of the target population who will recieve the
#' drug
#' @export
set_smc <- function(
  parameters,
  drug,
  start,
  end,
  frequency,
  min_age,
  max_age,
  coverage
  ) {
  parameters$smc <- TRUE
  parameters$smc_drug <- drug
  parameters$smc_start <- start
  parameters$smc_end <- end
  parameters$smc_frequency <- frequency
  parameters$smc_min_age <- min_age
  parameters$smc_max_age <- max_age
  parameters$smc_coverage <- coverage
  parameters
}

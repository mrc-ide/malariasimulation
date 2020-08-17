#' @title Parameterise an RTS,S strategy
#' @param parameters a list of parameters to modify
#' @param start the timestep to start
#' @param end the last timestep for the intervention
#' @param frequency the number of timsteps between rounds
#' @param min_ages for the target population, inclusive (in timesteps)
#' @param max_ages for the target population, inclusive (in timesteps)
#' @param boosters the timesteps (following the initial vaccination) at which booster vaccinations are administered
#' @param coverage the proportion of the target population who will recieve the
#' vaccine
#' @export
set_rtss <- function(
  parameters,
  start,
  end,
  frequency,
  min_ages,
  max_ages,
  boosters,
  coverage
  ) {
  parameters$rtss <- TRUE
  parameters$rtss_start <- start
  parameters$rtss_end <- end
  parameters$rtss_frequency <- frequency
  parameters$rtss_min_ages <- min_ages
  parameters$rtss_max_ages <- max_ages
  parameters$rtss_boosters <- boosters
  parameters$rtss_coverage <- coverage
  parameters
}

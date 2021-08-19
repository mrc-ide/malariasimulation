#' @title Parameterise clinical treatment
#'
#' @param parameters the model parameters
#' @param agegroups vector of agegroups (in timesteps)
#' @param timesteps vector of timesteps for each change in demography
#' @param birthrates vector of birthrates for each change in demography
#' @param deathrates matrix of deathrates per age group per timestep.
#' Rows are timesteps from the `timesteps` param. Columns are the age groups
#' from the `agegroups` param.
#' @export
set_demography <- function(
  parameters,
  agegroups,
  timesteps,
  birthrates,
  deathrates
  ) {
  stopifnot(all(agegroups > 0))
  stopifnot(all(timesteps > 0))
  stopifnot(all(deathrates > 0 & deathrates < 1))
  stopifnot(length(agegroups) == ncol(deathrates))
  stopifnot(length(timesteps) == nrow(deathrates))
  parameters$demography_agegroups <- agegroups
  parameters$demography_timesteps <- timesteps
  parameters$demography_birthrates <- birthrates
  parameters$demography_deathrates <- deathrates
  
  parameters
}

#' @title Parameterise clinical treatment
#'
#' @param parameters the model parameters
#' @param timesteps vector of timesteps for each change in population
#' @param population vector of new populations
set_population <- function(parameters, timesteps, population) {
  stopifnot(all(timesteps > 0))
  stopifnot(all(population > 0))
  parameters$population_timesteps <- timesteps
  parameters$human_population <- population

  parameters
}

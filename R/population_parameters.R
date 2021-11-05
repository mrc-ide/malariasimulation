#' @title Parameterise variable deathrates
#'
#' @param parameters the model parameters
#' @param agegroups vector of agegroups (in timesteps)
#' @param timesteps vector of timesteps for each change in demography
#' @param deathrates matrix of deathrates per age group per timestep.
#' Rows are timesteps from the `timesteps` param. Columns are the age groups
#' from the `agegroups` param.
#' @export
set_deathrates <- function(
  parameters,
  agegroups,
  timesteps,
  deathrates
  ) {
  stopifnot(all(agegroups > 0))
  stopifnot(all(timesteps > 0))
  stopifnot(all(deathrates > 0 & deathrates < 1))
  stopifnot(length(agegroups) == ncol(deathrates))
  stopifnot(length(timesteps) == nrow(deathrates))
  parameters$deathrate_agegroups <- agegroups
  parameters$deathrate_timesteps <- timesteps
  parameters$deathrates <- deathrates
  
  parameters
}

#' @title Parameterise variable birthrates
#'
#' @param parameters the model parameters
#' @param timesteps vector of timesteps for each change in birthrate
#' @param birthrates vector of birthrates for each timestep
set_birthrates <- function(parameters, timesteps, birthrates) {
  stopifnot(all(timesteps > 0))
  stopifnot(all(birthrates > 0))
  stopifnot(length(birthrates) == length(timesteps))
  parameters$birthrate_timesteps <- timesteps
  parameters$birthrates <- birthrates

  parameters
}

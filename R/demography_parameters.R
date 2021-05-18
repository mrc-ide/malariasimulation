#' @title Parameterise clinical treatment
#'
#' @param parameters the model parameters
#' @param agegroups vector of agegroups
#' @param timesteps vector of timesteps for each change in demography
#' @param birthrates vector of birthrates for each change in demography
#' @param deathrates matrix of deathrates per age group per timestep.
#' Rows are timesteps from the `timesteps` param. Columns are the age groups
#' from the `agegroups` param.
#' @export
set_demography <- function(parameters, agegroups, timesteps, birthrates, deathrates) {
  
  parameters$demography_agegroups <- agegroups # these should be in years
  parameters$demography_timesteps <- timesteps
  parameters$demography_birthrates <- birthrates
  parameters$demography_deathrates <- deathrates
  # add validation - age groups are +, timestamps are +, etc. error messages
  
  return(parameters)
  
}

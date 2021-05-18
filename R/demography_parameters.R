#' @title Parameterise clinical treatment
#'
#' @param parameters the model parameters
#' @param agegroups vector of agegroups
#' @param timesteps vector of timesteps for each change in demography
#' @param birthrates vector of birthrates for each change in demography
#' @param deathrates matrix of deathrates for each change in demography
#' @export
#' 
set_demography <- function(parameters, agegroups, timesteps, birthrates, deathrates) {
  
  parameters$demography_agegroups <- agegroups
  parameters$demography_timesteps <- timesteps
  parameters$demography_birthrates <- birthrates
  parameters$demography_deathrates <- deathrates
  # add validation - age groups are +, timestamps are +, etc. error messages
  
  return(parameters)
  
}
#' @title Get model parameters
#' @description
#' get_paramaters creates a list of parameters for use in the model. These
#' parameters are passed to process functions. These parameters are explained in
#' "The US President's Malaria Initiative, Plasmodium falciparum transmission
#' and mortality: A modelling study."
#'
#' NOTE: this function is likely to be extended to read in command line / config
#' file parameters
get_parameters <- function() {
  timestep_to_day <- 1
  human_population <- 100 * 1000
  parameters <- list(
    rd    = timestep_to_day / 5,
    ra    = timestep_to_day / 195,
    ru    = timestep_to_day / 110,
    rp    = timestep_to_day / 25,
    ft    = 1/2, # NOTE: set from sitefile
    rel   = 1 / (6.64 * timestep_to_day),
    rl    = 1 / (3.72 * timestep_to_day),
    rpl   = 1 / (.643 * timestep_to_day),
    mup   = .249,
    mum   = .249, #NOTE: set from sitefile
    beta  = 21.2,
    human_population = human_population,
    mosquito_limit   = 100 * human_population,
    timestep_to_day  = timestep_to_day
  )

  parameters
}

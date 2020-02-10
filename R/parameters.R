#' @title Get model parameters
#' @description
#' get_paramaters creates a named list of parameters for use in the model. These
#' parameters are passed to process functions. These parameters are explained in
#' "The US President's Malaria Initiative, Plasmodium falciparum transmission
#' and mortality: A modelling study."
#'
#' NOTE: this function is likely to be extended to read in command line / config
#' file parameters
#'
#' The parameters are defined as:
#' 
#' * days_per_timestep - the number of days to model per timestep
#' * rd - the rate at which humans move from state D to A
#' * ra - the rate at which humans move from state A to U
#' * ru - the rate at which humans move from state U to S
#' * rel - the rate at which mosquitos move from state E to L
#' * rl - the rate at which mosquitos move from state L to P
#' * rpl - the rate at which mosquitos move from state P to Sm
#' * mup - the rate at which pupal mosquitos die
#' * mum - the rate at which developed mosquitos die
#' * beta - the average number of eggs laid per female mosquito per day
#' * human population - the number of humans to model
#' * mosquito limit - the maximum number of mosquitos to allow for in the
#' * days_per_timestep - the number of days to model per timestep
get_parameters <- function() {
  days_per_timestep <- 1
  human_population <- 100 * 1000
  parameters <- list(
    rd    = days_per_timestep / 5,
    ra    = days_per_timestep / 195,
    ru    = days_per_timestep / 110,
    rel   = days_per_timestep / 6.64,
    rl    = days_per_timestep / 3.72,
    rpl   = days_per_timestep / .643,
    mup   = days_per_timestep * .249,
    mum   = days_per_timestep * .249, #NOTE: set from sitefile
    beta  = days_per_timestep * 21.2,
    human_population = human_population,
    mosquito_limit   = 100 * human_population,
    days_per_timestep  = days_per_timestep
  )

  parameters
}

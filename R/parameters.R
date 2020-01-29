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
  parameters <- list(
    rd    = 1/5,
    ra    = 1/195,
    ru    = 1/110,
    ft    = 1/2, # NOTE: set from sitefile
    rel   = 1 / (6.64 * timestep_to_day),
    rl    = 1 / (3.72 * timestep_to_day),
    rpl   = 1 / (.643 * timestep_to_day),
    mup   = .249,
    mum   = .249, #NOTE: set from sitefile
    timestep_to_day = timestep_to_day
  )

  parameters
}

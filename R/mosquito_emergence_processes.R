#' @title Mosquito emergence process
#' @description Move mosquitos from Unborn to Sm in line with the number of
#' pupals in the ODE models
#' 
#' @param mosquito the mosquito individual
#' @param odes a list of odes for each species
#' @param Unborn the unborn mosquito state
#' @param Sm the susceptable mosquito state
#' @param variety the variable representing the mosquito species
#' @param dpl the delay for pupal growth
#' @importFrom stats head
create_mosquito_emergence_process <- function(mosquito, odes, Unborn, Sm, variety, dpl) {
  function(api) {
    emergence_rate <- .5 * (1. - exp(-1./dpl))
    species <- c()
    for (v in seq_along(odes)) {
      n_pupals <- mosquito_model_get_states(odes[[v]])[[3]]
      species <- c(species, rep(v, n_pupals * emergence_rate))
    }
    if (length(species) > 0) {
      available <- api$get_state(mosquito, Unborn)
      if (length(available) < length(species)) {
        stop(paste0(
          'Run out of mosquitos. Short by ',
          length(species) - length(available),
          '. Please increase mosquito_limit'
        ))
      }
      target <- head(available, length(species))
      api$queue_state_update(mosquito, Sm, target)
      api$queue_variable_update(mosquito, variety, species, target)
    }
  }
}

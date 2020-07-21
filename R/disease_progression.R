#' @title Modelling the progression of the human disease
#' @description Randomly moves individuals towards the later stages of disease
#' and updates their infectivity
#'
#' @param human the handle for the human individuals
#' @param from_state the source disease state
#' @param to_state the destination disease state
#' @param rate the rate at which to move humans
#' @param infectivity the handle for the infectivity variable
#' @param new_infectivity the new infectivity of the progressed individuals
create_progression_process <- function(
 human,
 from_state,
 to_state,
 rate,
 infectivity,
 new_infectivity
 ) {
  function(api) {
    source_humans <- api$get_state(human, from_state)
    to_move <- source_humans[bernoulli(length(source_humans), rate)]
    if (length(to_move) > 0) {
      api$queue_state_update(human, to_state, to_move)
      api$queue_variable_update(human, infectivity, new_infectivity, to_move)
    }
  }
}

#' @title Modelling the progression to asymptomatic disease
#' @description Randomly moves individuals to asymptomatic disease and
#' calculates the infectivity for their age and immunity
#'
#' @param human the handle for the human individuals
#' @param states the available human states
#' @param variables the available human variables
#' @param rate the rate of progression to asymptomatic infection
create_asymptomatic_progression_process <- function(
  human,
  states,
  variables,
  rate
  ) {
  function(api) {
    source_humans <- api$get_state(human, states$D)
    to_move <- source_humans[bernoulli(length(source_humans), rate)]
    if (length(to_move) > 0) {
      api$queue_state_update(human, states$A, to_move)
      new_infectivity <- asymptomatic_infectivity(
        get_age(
          api$get_variable(human, variables$birth, to_move),
          api$get_timestep()
        ),
        api$get_variable(human, variables$id, to_move),
        api$get_parameters()
      )
      api$queue_variable_update(
        human,
        variables$infectivity,
        new_infectivity,
        to_move
      )
    }
  }
}

#' @title Update the state of an individual as infection events occur
#' @description Randomly moves individuals towards the later stages of disease
#' and updates their infectivity
#'
#' @param state the human state variable
#' @param to_state the destination disease state
#' @param infectivity the handle for the infectivity variable
#' @param new_infectivity the new infectivity of the progressed individuals
#' @noRd
create_infection_update_listener <- function(
 state,
 to_state,
 infectivity,
 new_infectivity
 ) {
  function(timestep, to_move) {
    state$queue_update(to_state, to_move)
    infectivity$queue_update(new_infectivity, to_move)
  }
}

create_progression_process <- function(event, state, from_state, rate) {
  function(timestep) {
    event$schedule(state$get_index_of(from_state)$sample(1/rate), 0)
  }
}

create_rate_listener <- function(from_state, to_state, renderer) {
  function(timestep, target) {
    renderer$render(
      paste0('rate_', from_state, '_', to_state),
      target$size(),
      timestep
    )
  }
}

#' @title Modelling the progression to asymptomatic disease
#' @description Randomly moves individuals to asymptomatic disease and
#' calculates the infectivity for their age and immunity
#'
#' @param variables the available human variables
#' @param parameters model parameters
#' @noRd
create_asymptomatic_update_listener <- function(variables, parameters) {
  function(timestep, to_move) {
    if (to_move$size() > 0) {
      variables$state$queue_update('A', to_move)
      new_infectivity <- asymptomatic_infectivity(
        get_age(
          variables$birth$get_values(to_move),
          timestep
        ),
        variables$id$get_values(to_move),
        parameters
      )
      variables$infectivity$queue_update(
        new_infectivity,
        to_move
      )
    }
  }
}

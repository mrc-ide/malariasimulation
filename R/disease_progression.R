#' @title Update the state of an individual as infection events occur
#' @description Randomly moves individuals towards the later stages of disease
#' and updates their infectivity
#'
#' @param state the human state variable
#' @param to_state the destination disease state
#' @param infectivity the handle for the infectivity variable
#' @param new_infectivity the new infectivity of the progressed individuals
#' @noRd
update_infection <- function(
  state,
  to_state,
  infectivity,
  new_infectivity,
  to_move
  ) {
  state$queue_update(to_state, to_move)
  infectivity$queue_update(new_infectivity, to_move)
}

create_progression_process <- function(
  state,
  from_state,
  to_state,
  rate,
  infectivity,
  new_infectivity
  ) {
  function(timestep) {
    to_move <- state$get_index_of(from_state)$sample(rate)
    update_infection(
      state,
      to_state,
      infectivity,
      new_infectivity,
      to_move
    )
  }
}

create_asymptomatic_progression_process <- function(
  state,
  rate,
  variables,
  parameters
  ) {
  function(timestep) {
    to_move <- state$get_index_of('D')$sample(rate)
    update_to_asymptomatic_infection(
      variables,
      parameters,
      timestep,
      to_move
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
update_to_asymptomatic_infection <- function(
  variables,
  parameters,
  timestep,
  to_move
  ) {
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

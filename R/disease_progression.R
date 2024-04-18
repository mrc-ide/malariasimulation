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
    
    # Retrieve the indices of all individuals in the to_move state:
    index <- state$get_index_of(from_state)
    
    # If the length of rate is greater than 1 (when it's a variable):
    if (inherits(rate, "DoubleVariable")) {
      rate <- rate$get_values(index)
    }
    
    # Sample the individuals to be moved into a new Bitset using the transition rate(s):
    to_move <- index$sample(1/rate)
    
    # Update the infection status of those individuals who are moving:
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
    to_move <- state$get_index_of('D')$sample(1/rate)
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

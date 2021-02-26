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

#' @title Schedule progression of human disease at the start of the simulation
#'
#' @param event the event to schedule
#' @param state the human state variable
#' @param from_state the state this event applies to
#' @param rate the average time spent in this state
#' @noRd
initialise_progression <- function(event, state, from_state, rate) {
  target <- state$get_index_of(from_state)
  event$schedule(
    target,
    log_uniform(target$size(), rate)
  )
}

#' @title Modelling the progression of the human disease
#' @description schedules follow up infection events with the log uniform dist.
#'
#' @param event the event to schedule
#' @param rate the average time spent in this state
#' @noRd
create_progression_listener <- function(event, rate) {
  function(timestep, target) {
    event$schedule(target, log_uniform(target$size(), rate))
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

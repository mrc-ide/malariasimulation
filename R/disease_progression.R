#' @title Update the state of an individual as infection events occur
#' @description Randomly moves individuals towards the later stages of disease
#' and updates their infectivity
#'
#' @param human the handle for the human individuals
#' @param to_state the destination disease state
#' @param infectivity the handle for the infectivity variable
#' @param new_infectivity the new infectivity of the progressed individuals
create_infection_update_listener <- function(
 human,
 to_state,
 infectivity,
 new_infectivity
 ) {
  function(api, to_move) {
    api$queue_state_update(human, to_state, to_move)
    api$queue_variable_update(human, infectivity, new_infectivity, to_move)
  }
}

#' @title Schedule progression of human disease at the start of the simulation
#'
#' @param api the simulation api
#' @param event the event to schedule
#' @param human the human handle
#' @param from_state the state this event applies to
#' @param rate the average time spent in this state
initialise_progression <- function(api, event, human, from_state, rate) {
  target <- api$get_state(human, from_state)
  api$schedule(event, target, log_uniform(length(target), rate))
}

#' @title Modelling the progression of the human disease
#' @description schedules follow up infection events with the log uniform dist.
#'
#' @param event the event to schedule
#' @param rate the average time spent in this state
create_progression_listener <- function(event, rate) {
  function(api, target) {
    api$schedule(event, target, log_uniform(length(target), rate))
  }
}

#' @title Modelling the progression to asymptomatic disease
#' @description Randomly moves individuals to asymptomatic disease and
#' calculates the infectivity for their age and immunity
#'
#' @param human the handle for the human individuals
#' @param states the available human states
#' @param variables the available human variables
create_asymptomatic_update_listener <- function(
  human,
  states,
  variables
  ) {
  function(api, to_move) {
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

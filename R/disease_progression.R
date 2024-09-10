#' @title Calculate recovery rates
#' @description Calculates recovery rates for each individual in the population
#' for storage in competing hazards object and subsequent resolution
#'
#' @param variables the available human variables
#' @param recovery_outcome competing hazards object for recovery rates
#' @noRd
create_recovery_rates_process <- function(
  variables,
  recovery_outcome
) {
  function(timestep){
    target <- variables$state$get_index_of("S")$not()
    recovery_outcome$set_rates(
      target,
      variables$recovery_rates$get_values(target))
  }
}


#' @title Disease progression outcomes (recovery)
#' @description Following resolution of competing hazards, update state and
#' infectivity of sampled individuals
#'
#' @param timestep the current timestep
#' @param target the sampled recovering individuals
#' @param variables the available human variables
#' @param parameters model parameters
#' @param renderer competing hazards object for recovery rates
#' @noRd
recovery_outcome_process <- function(
    timestep,
    target,
    variables,
    parameters,
    renderer
){
  
  update_to_asymptomatic_infection(
    variables,
    parameters,
    timestep,
    variables$state$get_index_of("D")$and(target)
  )
  
  update_infection(
    variables$state,
    "U",
    variables$infectivity,
    parameters$cu,
    variables$recovery_rates,
    1/parameters$du,
    variables$state$get_index_of("A")$and(target)
  )
  
  update_infection(
    variables$state,
    "S",
    variables$infectivity,
    0,
    variables$recovery_rates,
    0,
    variables$state$get_index_of(c("U","Tr"))$and(target)
  )
  
}

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
    recovery_rates,
    new_recovery_rate,
    to_move
) {
  state$queue_update(to_state, to_move)
  infectivity$queue_update(new_infectivity, to_move)
  recovery_rates$queue_update(new_recovery_rate, to_move)
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
    variables$recovery_rates$queue_update(
      1/parameters$da,
      to_move
    )
  }
}

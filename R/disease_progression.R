#' @title Calculate recovery rates
#' @description Calculates recovery rates for each individual in the population
#' for storage in competing hazards object and future resolution
#'
#' @param variables the available human variables
#' @param parameters model parameters
#' @param recovery_outcome competing hazards object for recovery rates
#' @noRd
calculate_recovery_rates <- function(variables, parameters, recovery_outcome, dt_input){
  
  # Get correct input for dt depending on antimalarial resistance
  if(isFALSE(parameters$antimalarial_resistance) & length(dt_input) == 1 & is.numeric(dt_input)){
    dt_v <- dt_input
  } else if (isTRUE(parameters$antimalarial_resistance)) {
    dt_v <- dt_input$get_values(variables$state$get_index_of("Tr"))
  }

  recovery_rates <- numeric(length = parameters$human_population)
  recovery_rates[variables$state$get_index_of("D")$to_vector()] <- 1/parameters$dd
  recovery_rates[variables$state$get_index_of("A")$to_vector()] <- 1/parameters$da
  recovery_rates[variables$state$get_index_of("U")$to_vector()] <- 1/parameters$du
  recovery_rates[variables$state$get_index_of("Tr")$to_vector()] <- 1/dt_v
  recovery_outcome$set_rates(recovery_rates)
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
recovery_process_resolved_hazard <- function(
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
    variables$state$get_index_of("A")$and(target)
  )
  
  update_infection(
    variables$state,
    "S",
    variables$infectivity,
    0,
    variables$state$get_index_of("U")$and(target)
  )
  
  update_infection(
    variables$state,
    "S",
    variables$infectivity,
    0,
    variables$state$get_index_of("Tr")$and(target)
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
  to_move
  ) {
  state$queue_update(to_state, to_move)
  infectivity$queue_update(new_infectivity, to_move)
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

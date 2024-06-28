#' @title Calculate disease progression rates
#' @description Calculates disease progression rates for each individual in the population
#' for storage in competing hazards object and subsequent resolution
#'
#' @param variables the available human variables
#' @param progression_outcome competing hazards object for disease progression rates
#' @noRd
create_progression_rates_process <- function(
  variables,
  progression_outcome
) {
  function(timestep){
    target <- variables$state$get_index_of("S")$not()
    progression_outcome$set_rates(
      target,
      variables$progression_rates$get_values(target))
  }
}


#' @title Disease progression outcomes
#' @description Following resolution of competing hazards, update state and
#' infectivity of sampled individuals
#'
#' @param timestep the current timestep
#' @param target the sampled progressing individuals
#' @param variables the available human variables
#' @param parameters model parameters
#' @param renderer competing hazards object for disease progression rates
#' @noRd
progression_outcome_process <- function(
    timestep,
    target,
    variables,
    parameters,
    renderer
){
  
  if(parameters$parasite == "falciparum"){
    # p.f has immunity-determined asymptomatic infectivity
    update_to_asymptomatic_infection(
      variables,
      parameters,
      timestep,
      variables$state$get_index_of("D")$and(target)
    )
  } else if (parameters$parasite == "vivax"){
    # p.v has constant asymptomatic infectivity
    update_infection(
      variables$state,
      "A",
      variables$infectivity,
      parameters$ca,
      variables$progression_rates,
      1/parameters$da,
      variables$state$get_index_of("D")$and(target)
    )
  }
  
  update_infection(
    variables$state,
    "U",
    variables$infectivity,
    parameters$cu,
    variables$progression_rates,
    1/parameters$du,
    variables$state$get_index_of("A")$and(target)
  )
  
  update_infection(
    variables$state,
    "S",
    variables$infectivity,
    0,
    variables$progression_rates,
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
#' @param progression_rates the handle for the progression_rates variable
#' @param new_progression the new disease progression rate of the progressed individuals
#' @noRd
update_infection <- function(
    state,
    to_state,
    infectivity,
    new_infectivity,
    progression_rates,
    new_progression_rate,
    to_move
) {
  state$queue_update(to_state, to_move)
  infectivity$queue_update(new_infectivity, to_move)
  progression_rates$queue_update(new_progression_rate, to_move)
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
    variables$progression_rates$queue_update(
      1/parameters$da,
      to_move
    )
  }
}



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
    api$queue_state_update(human, to_state, to_move)
    api$queue_variable_update(human, infectivity, new_infectivity, to_move)
  }
}

create_asymptomatic_progression_process <- function(
  human,
  states,
  variables,
  rate
  ) {
  source_humans <- api$get_state(human, states$D)
  to_move <- source_humans[bernoulli(length(source_humans), rate)]
  api$queue_state_update(human, states$A, to_move)
  new_infectivity <- asymptomatic_infectivity(
    get_age(api$get_variable(human, variables$birth, to_move), api$get_timestep()),
    api$get_variable(human, variables$id, to_move),
    api$get_parameters()
  )
  api$queue_variable_update(human, infectivity, new_infectivity, to_move)
}

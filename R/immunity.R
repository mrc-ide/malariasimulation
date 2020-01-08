maternal_immunity_process <- function(simulation_frame, timestep, parameters) {
  i <- simulation_frame$get_variable(human, icm)
  VariableUpdate$new(human, icm, i - parameters$rm * i)
}

preerythoctic_immunity_process <- function(simulation_frame, timestep, parameters) {
  i <- simulation_frame$get_variable(human, ib)
  VariableUpdate$new(
    human,
    ib,
    immunity_decay(
      i,
      simulation_frame$get_variable(human, last_bitten),
      timestep,
      parameters$rb,
      parameters$ub
    )
  )
}

acquired_immunity_process <- function(simulation_frame, timestep, parameters) {
  i <- simulation_frame$get_variable(human, ica)
  VariableUpdate$new(
    human,
    ica,
    immunity_decay(
      i,
      simulation_frame$get_variable(human, last_infected),
      timestep,
      parameters$rc,
      parameters$uc
    )
  )
}

immunity_decay <- function(level, last_timestep, timestep, rate, delay) {
  boost <- (timestep - last_timestep) == delay
  level[boost] <- level[boost] + 1
  level[!boost] <- level[!boost] - rate * level[!boost]
  level
}

egg_laying_process <- function(simulation_frame, timestep, parameters) {
  m <- simulation_frame$get_state(mosquito, Sm, Im)
  unborn <- simulation_frame$get_state(mosquito, Unborn)
  if (length(m) > 0) {
    n_eggs <- parameters$beta * length(m)
    if (n_eggs > length(unborn)) {
      stop('Run out of mosquitos')
    }
    return(individual::StateUpdate$new(mosquito, E, unborn[seq_len(n_eggs)]))
  }
}

larval_death_process <- function(simulation_frame, timestep, parameters) {
  early_larval <- simulation_frame$get_state(mosquito, E)
  late_larval <- simulation_frame$get_state(mosquito, L)
  n <- length(early_larval) + length(late_larval)
  k <- carrying_capacity(timestep, parameters)
  early_regulation <- 1 + n / k
  late_regulation <- 1 + parameters$gamma * n / k
  early_larval_deaths <- early_larval[
    runif(length(early_larval), 0, 1) < parameters$me * early_regulation
  ]
  late_larval_deaths <- late_larval[
    runif(length(late_larval), 0, 1) < parameters$ml * late_regulation
  ]
  list(
    individual::StateUpdate$new(mosquito, Unborn, early_larval_deaths),
    individual::StateUpdate$new(mosquito, Unborn, late_larval_deaths)
  )
}

carrying_capacity <- function(timestep, parameters) {
  rainfall <- parameters$g0 + sum(vnapply(seq_len(3), function(i) {
    parameters[[
      'g_' + i
    ]] * cos(2 * pi * timestep * i) + parameters[[
      'h_' + i
    ]] * sin(2 * pi * timestep * i)
  }))
  parameters$K0 * rainfall / parameters$R_bar
}

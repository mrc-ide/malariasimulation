#' @title Mosquito births
#' @description
#' This is the process for mosquito birth, it defines how many new early stage
#' larvae are created on each timestep.
#' @param mosquito, the mosquito individual
#' @param Sm, the susceptable mosquito state
#' @param Im, the infected mosquito state
#' @param Unborn, the unborn mosquito state
#' @param E, the early stage larval state
create_egg_laying_process <- function(mosquito, Sm, Im, Unborn, E) {
  function(simulation_frame, timestep, parameters) {
    m <- simulation_frame$get_state(mosquito, Sm, Im)
    unborn <- simulation_frame$get_state(mosquito, Unborn)
    if (length(m) > 0) {
      n_eggs <- parameters$beta * length(m)
      if (n_eggs > length(unborn)) {
        stop('Run out of mosquitos')
      }
      if (n_eggs >= 1) {
        return(individual::StateUpdate$new(mosquito, E, unborn[seq_len(n_eggs)]))
      }
    }
  }
}

#' @title Larval deaths
#' @description
#' This process defines how many early and late stage larvae die due to
#' seasonal carrying capacity.
#' @param mosquito, the mosquito individual
#' @param E, the early stage larval state
#' @param L, the late stage larval state
#' @param Unborn, the unborn mosquito state
create_larval_death_process <- function(mosquito, E, L, Unborn) {
  function(simulation_frame, timestep, parameters) {
    early_larval <- simulation_frame$get_state(mosquito, E)
    late_larval <- simulation_frame$get_state(mosquito, L)
    n <- length(early_larval) + length(late_larval)
    k <- carrying_capacity(timestep, parameters)
    early_regulation <- 1 + n / k
    late_regulation <- 1 + parameters$gamma * n / k
    early_larval_deaths <- early_larval[
      bernoulli(length(early_larval), parameters$me * early_regulation)
    ]
    late_larval_deaths <- late_larval[
      bernoulli(length(late_larval), parameters$ml * late_regulation)
    ]
    individual::StateUpdate$new(
      mosquito,
      Unborn,
      c(early_larval_deaths, late_larval_deaths)
    )
  }
}

carrying_capacity <- function(timestep, parameters) {
  r <- rainfall(
    timestep,
    parameters$days_per_timestep,
    parameters$g0,
    c(parameters$g1, parameters$g2, parameters$g3),
    c(parameters$h1, parameters$h2, parameters$h3)
  )
  parameters$K0 * r / parameters$R_bar
}

rainfall <- function(t, days_per_timestep, g0, g, h) {
  g0 + sum(vnapply(seq_len(3), function(i) {
    g[i] * cos(
      2 * pi * t * 365 / days_per_timestep * i
    ) + h[i] * sin(
      2 * pi * t * 365 / days_per_timestep * i
    )
  }))
}

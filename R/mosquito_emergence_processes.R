#' @title Larval deaths
#' @description
#' This process defines how many early and late stage larvae die due to
#' seasonal carrying capacity.
#' @param mosquito the mosquito individual
#' @param E the early stage larval state
#' @param L the late stage larval state
#' @param Unborn the unborn mosquito state
#' @param events a list of events in the model
create_larval_death_process <- function(mosquito, E, L, Unborn, events) {
  function(api) {
    timestep <- api$get_timestep()
    parameters <- api$get_parameters()
    early_larval <- api$get_state(mosquito, E)
    late_larval <- api$get_state(mosquito, L)
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
    api$clear_schedule(events$larval_growth, early_larval_deaths)
    api$clear_schedule(events$pupal_development, late_larval_deaths)
    api$queue_state_update(
      mosquito,
      Unborn,
      c(early_larval_deaths, late_larval_deaths)
    )
  }
}

create_pupal_death_process <- function(mosquito, P, Unborn, rate, events) {
  function (api) {
    source_individuals <- api$get_state(mosquito, P)
    target_individuals <- source_individuals[
      bernoulli(length(source_individuals), rate)
    ]
    api$clear_schedule(events$susceptable_development, target_individuals)
    api$queue_state_update(mosquito, Unborn, target_individuals)
  }
}

create_mosquito_death_process <- function(mosquito, states, rate, events) {
  function (api) {
    source_individuals <- api$get_state(
      mosquito,
      states$Sm,
      states$Pm,
      states$Im
    )
    target_individuals <- source_individuals[
      bernoulli(length(source_individuals), rate)
    ]
    api$clear_schedule(events$mosquito_infection, target_individuals)
    api$queue_state_update(mosquito, states$Unborn, target_individuals)
  }
}

carrying_capacity <- function(timestep, parameters) {
  if (parameters$model_seasonality) {
    r <- rainfall(
      timestep,
      parameters$days_per_timestep,
      parameters$g0,
      c(parameters$g1, parameters$g2, parameters$g3),
      c(parameters$h1, parameters$h2, parameters$h3)
    )
    return(parameters$K0 * r / parameters$R_bar)
  }
  parameters$K0
}

rainfall <- function(t, days_per_timestep, g0, g, h) {
  g0 + sum(vnapply(seq_len(3), function(i) {
    g[i] * cos(
      2 * pi * t * days_per_timestep / 365 * i
    ) + h[i] * sin(
      2 * pi * t * days_per_timestep / 365 * i
    )
  }))
}

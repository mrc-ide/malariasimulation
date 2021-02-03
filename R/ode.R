ODE_INDICES <- c(E = 1, L = 2, P = 3)

parameterise_ode <- function(parameters) {
  lapply(
    parameters$species_proportions,
    function(p) {
      m <- p * parameters$total_M
      create_mosquito_model(
        initial_mosquito_counts(parameters, 0, m)[ODE_INDICES],
        parameters$beta,
        parameters$del,
        parameters$me,
        calculate_carrying_capacity(parameters),
        parameters$gamma,
        parameters$dl,
        parameters$ml,
        parameters$dpl,
        parameters$mup,
        m,
        parameters$model_seasonality,
        parameters$days_per_timestep,
        parameters$g0,
        parameters$g,
        parameters$h,
        calculate_R_bar(parameters)
      )
    }
  )
}

create_ode_rendering_process <- function(renderer, odes) {
  function(timestep) {
    counts <- rep(0, length(ODE_INDICES))
    for (ode in odes) {
      row <- mosquito_model_get_states(ode)
      counts <- counts + row
    }
    for (i in seq_along(ODE_INDICES)) {
      renderer$render(
        paste0('mosquito_', names(ODE_INDICES)[[i]], '_count'),
        counts[[i]],
        timestep
      )
    }
  }
}

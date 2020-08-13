parameterise_ode <- function(parameters) {
  lapply(
    parameters$variety_proportions,
    function(p) {
      m <- p * parameters$total_M
      create_mosquito_model(
        initial_mosquito_counts(parameters, 0, m)[seq(3)],
        parameters$beta,
        parameters$del,
        parameters$me,
        calculate_carrying_capacity(parameters),
        parameters$gamma,
        parameters$dl,
        parameters$ml,
        parameters$dpl,
        parameters$mup,
        m
      )
    }
  )
}

create_ode_rendering_process <- function(odes) {
  mosquito_states <- c('E', 'L', 'P')
  function(api) {
    counts <- rep(0, length(mosquito_states))
    for (ode in odes) {
      row <- mosquito_model_get_states(ode)
      counts <- counts + row
    }
    for (i in seq_along(mosquito_states)) {
      api$render(
        paste0('mosquito_', mosquito_states[[i]], '_count'),
        counts[[i]]
      )
    }
  }
}

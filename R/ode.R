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


#' @title Step mosquito ODE
#' @description collects summarises the human state, sends it to the vector ode
#' and makes a step
#'
#' @param odes the models to step, one for each species
#' @param mosquito the mosquito individual handle
#' @param states a list of all of the model states
#' @param variables a list of all of the model variables
create_ode_stepping_process <- function(
  odes,
  mosquito,
  states,
  variables
  ) {
  function(api) {
    variety <- api$get_variable(mosquito, variables$mosquito_variety)
    adult_mosquitoes <- api$get_state(
      mosquito,
      states$Sm,
      states$Pm,
      states$Im
    )
    for (species in seq_along(api$get_parameters()$blood_meal_rates)) {
      total_M <- length(intersect(which(variety == species), adult_mosquitoes))
      api$render(paste0('total_M_', species), total_M)
      mosquito_model_step(odes[[species]], total_M)
    }
  }
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

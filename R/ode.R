parameterise_ode <- function(parameters, foim = 0.) {
  lapply(
    parameters$variety_proportions,
    function(p) {
      m <- p * parameters$human_population * parameters$density
      create_mosquito_model(
        initial_mosquito_counts(parameters, foim, m),
        parameters$beta,
        parameters$del,
        parameters$me,
        parameters$K0,
        parameters$gamma,
        parameters$dl,
        parameters$ml,
        parameters$dpl,
        parameters$mup,
        foim,
        parameters$mum,
        parameters$dem
      )
    }
  )
}


#' @title Step mosquito ODE
#' @description collects summarises the human state, sends it to the vector ode
#' and makes a step
#'
#' @param ode the models to step, one for each species
#' @param human the human individual
#' @param states a list of all of the model states
#' @param variables a list of all of the model variables
create_ode_stepping_process <- function(
  odes,
  human,
  states,
  variables
  ) {
  function(api) {
    parameters <- api$get_parameters()
    age <- get_age(api$get_variable(human, variables$birth), api$get_timestep())
    xi  <- api$get_variable(human, variables$xi)
    a_subset <- api$get_state(human, states$A)
    d_subset <- api$get_state(human, states$D)
    u_subset <- api$get_state(human, states$U)

    a_infectivity <- asymptomatic_infectivity(
      age[a_subset],
      api$get_variable(human, variables$id, a_subset),
      parameters
    )

    age_subset <- c(age[d_subset], age[a_subset], age[u_subset])
    xi_subset  <- c(xi[d_subset], xi[a_subset], xi[u_subset])

    infectivity<- c(
      rep(parameters$cd, length(d_subset)),
      a_infectivity,
      rep(parameters$cu, length(u_subset))
    )

    lambda <- mosquito_force_of_infection(
      seq_along(parameters$blood_meal_rates),
      age_subset,
      xi_subset,
      infectivity,
      parameters
    )

    for (species in seq_along(lambda)) {
      api$render(paste0('FOIM_', species), lambda[[species]])
      mosquito_model_step(odes[[species]], lambda[[species]])
    }
  }
}

create_ode_rendering_process <- function(odes) {
  mosquito_states <- c(
    'E',
    'L',
    'P',
    'Sm',
    'Pm',
    'Im'
  )
  function(api) {
    counts <- rep(0, length(mosquito_states))
    for (ode in odes) {
      row <- mosquito_model_get_states(ode)
      counts <- counts + row
    }
    for (i in seq_along(mosquito_states)) {
      api$render(
        paste0('mosquito_', mosquito_states[[i]], '_counts'),
        counts[[i]]
      )
    }
  }
}

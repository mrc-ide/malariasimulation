ODE_INDICES <- c(E = 1, L = 2, P = 3)
ADULT_ODE_INDICES <- c(Sm = 4, Pm = 5, Im = 6)

parameterise_mosquito_models <- function(parameters) {
  lapply(
    parameters$species_proportions,
    function(p) {
      m <- p * parameters$total_M
      growth_model <- create_mosquito_model(
        parameters$beta,
        parameters$del,
        parameters$me,
        p * calculate_carrying_capacity(parameters),
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

      if (!parameters$hybrid_mosquitoes) {
        susceptible <- initial_mosquito_counts(parameters, 0, m)[ADULT_ODE_INDICES['Sm']]
        return(
          create_adult_mosquito_model(
            growth_model,
            parameters$mum,
            parameters$dem,
            susceptible
          )
        )
      }
      growth_model
    }
  )
}

parameterise_solvers <- function(models, parameters) {
  lapply(
    seq_along(models),
    function(i) {
      m <- parameters$species_proportions[[i]] * parameters$total_M
      if (!parameters$hybrid_mosquitoes) {
        return(
          create_adult_solver(
            models[[i]], 
            initial_mosquito_counts(parameters, 0, m)
          )
        )
      }
      create_solver(
        models[[i]], 
        initial_mosquito_counts(parameters, 0, m)[ODE_INDICES]
      )
    }
  )
}

create_ode_rendering_process <- function(renderer, solvers) {
  if (parameters$hybrid_mosquitoes) {
    indices <- c(ODE_INDICES, ADULT_ODE_INDICES)
  } else {
    indices <- ODE_INDICES
  }

  function(timestep) {
    counts <- rep(0, length(indices))
    for (ode in odes) {
      row <- solver_get_states(ode)
      counts <- counts + row
    }
    for (i in seq_along(indices)) {
      renderer$render(
        paste0(names(indices)[[i]], '_count'),
        counts[[i]],
        timestep
      )
    }
  }
}

#' @title Step mosquito solver
#' @description calculates total_M per species and updates the vector ode
#'
#' @param models for each species
#' @param solvers for each species
#' @param state the mosquito state variable
#' @param species the mosquito species variable
#' @param species_names the names of the mosquito species
#' @param renderer the model renderer
#' @noRd
create_solver_stepping_process <- function(
  modles,
  solvers,
  state,
  species,
  species_names,
  renderer
  ) {
  function(timestep) {
    adult <- state$get_index_of("NonExistent")$not()
    for (s_i in seq_along(species_names)) {
      total_M <- species$get_index_of(species_names[[s_i]])$and(adult)$size()
      renderer$render(paste0('total_M_', s_i), total_M, timestep)
      mosquito_model_update(models[[s_i]], total_M)
      solver_step(solvers[[s_i]])
    }
    renderer$render('total_M', adult$size(), timestep)
  }
}


#' @title Step adult mosquito solver
#' @description steps the full vector ODE (not just growth)
#'
#' @param solvers the models to step, one for each species
#' @param renderer the model renderer
#' @noRd
create_adult_stepping_process <- function(
  solvers,
  renderer
  ) {
  function(timestep) {
    total_M <- 0
    for (s_i in seq_along(solvers)) {
      solver_step(solvers[[s_i]])
      species_total_M <- sum(
        solver_get_states(solvers[[s_i]])[ADULT_ODE_INDICES]
      )
      renderer$render(paste0('total_M_', s_i), species_total_M, timestep)
      total_M <- total_M + species_total_M
    }
    renderer$render('total_M', total_M, timestep)
  }
}

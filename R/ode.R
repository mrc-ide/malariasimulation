ODE_INDICES <- c(E = 1, L = 2, P = 3)
ADULT_ODE_INDICES <- c(Sm = 4, Pm = 5, Im = 6)

parameterise_mosquito_models <- function(parameters) {
  lapply(
    seq_along(parameters$species),
    function(i) {
      p <- parameters$species_proportions[[i]]
      m <- p * parameters$total_M
      growth_model <- create_mosquito_model(
        parameters$beta,
        parameters$del,
        parameters$me,
        calculate_carrying_capacity(parameters, m, i),
        parameters$gamma,
        parameters$dl,
        parameters$ml,
        parameters$dpl,
        parameters$mup,
        m,
        parameters$model_seasonality,
        parameters$g0,
        parameters$g,
        parameters$h,
        calculate_R_bar(parameters)
      )

      if (!parameters$hybrid_mosquitoes) {
        susceptible <- initial_mosquito_counts(
          parameters,
          i,
          parameters$init_foim,
          m
        )[ADULT_ODE_INDICES['Sm']]
        return(
          create_adult_mosquito_model(
            growth_model,
            parameters$mum[[i]],
            parameters$dem,
            susceptible * parameters$init_foim
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
      init <- initial_mosquito_counts(parameters, i, parameters$init_foim, m)
      if (!parameters$hybrid_mosquitoes) {
        return(
          create_adult_solver(models[[i]], init)
        )
      }
      create_solver(models[[i]], init[ODE_INDICES])
    }
  )
}

create_ode_rendering_process <- function(renderer, solvers, parameters) {
  if (parameters$hybrid_mosquitoes) {
    indices <- ODE_INDICES
  } else {
    indices <- c(ODE_INDICES, ADULT_ODE_INDICES)
  }

  function(timestep) {
    counts <- rep(0, length(indices))
    for (i in seq_along(solvers)) {
      row <- solver_get_states(solvers[[i]])
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
#' @param solvers for each species
#' @noRd
create_solver_stepping_process <- function(solvers) {
  function(timestep) {
    for (solver in solvers) {
      solver_step(solver)
    }
  }
}

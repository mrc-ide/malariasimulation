ODE_INDICES <- c(E = 1, L = 2, P = 3)
ADULT_ODE_INDICES <- c(Sm = 4, Pm = 5, Im = 6)

parameterise_mosquito_models <- function(parameters, timesteps) {
  
  lapply(
    seq_along(parameters$species),
    function(i) {
      p <- parameters$species_proportions[[i]]
      m <- p * parameters$total_M
      # Baseline carrying capacity
      k0 <- calculate_carrying_capacity(parameters, m, i)
      # Create the carrying capacity object
      k_timeseries <- create_timeseries(size = length(parameters$carrying_capacity_timesteps), k0)
      if(parameters$carrying_capacity){
        for(j in 1:length(parameters$carrying_capacity_timesteps)){
          timeseries_push(
            k_timeseries,
            parameters$carrying_capacity_scalers[j,i] * k0,
            parameters$carrying_capacity_timesteps[j]
          )
        }
      }
      growth_model <- create_aquatic_mosquito_model(
        parameters$beta,
        parameters$del,
        parameters$me,
        k_timeseries,
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
        calculate_R_bar(parameters),
        parameters$mum[[i]],
        parameters$blood_meal_rates[[i]],
        parameters$rainfall_floor
      )
      
      if (!parameters$individual_mosquitoes) {
        susceptible <- initial_mosquito_counts(
          parameters,
          i,
          parameters$init_foim,
          m
        )[ADULT_ODE_INDICES['Sm']]
        return(
          AdultMosquitoModel$new(create_adult_mosquito_model(
            growth_model,
            parameters$mum[[i]],
            parameters$dem,
            susceptible * parameters$init_foim,
            parameters$init_foim
          ))
        )
      }
      AquaticMosquitoModel$new(growth_model)
    }
  )
}

parameterise_solvers <- function(models, parameters) {
  lapply(
    seq_along(models),
    function(i) {
      m <- parameters$species_proportions[[i]] * parameters$total_M
      init <- initial_mosquito_counts(parameters, i, parameters$init_foim, m)
      if (!parameters$individual_mosquitoes) {
        return(
          Solver$new(create_adult_solver(
            models[[i]]$.model,
            init,
            parameters$r_tol,
            parameters$a_tol,
            parameters$ode_max_steps
          ))
        )
      }
      Solver$new(create_aquatic_solver(
        models[[i]]$.model,
        init[ODE_INDICES],
        parameters$r_tol,
        parameters$a_tol,
        parameters$ode_max_steps
      ))
    }
  )
}

create_compartmental_rendering_process <- function(renderer, solvers, parameters) {
  if (parameters$individual_mosquitoes) {
    indices <- ODE_INDICES
  } else {
    indices <- c(ODE_INDICES, ADULT_ODE_INDICES)
  }
  
  function(timestep) {
    counts <- rep(0, length(indices))
    for (s_i in seq_along(solvers)) {
      if (parameters$species_proportions[[s_i]] > 0) {
        row <- solvers[[s_i]]$get_states()
      } else {
        row <- rep(0, length(indices))
      }
      for (i in seq_along(indices)) {
        renderer$render(
          paste0(names(indices)[[i]], '_', parameters$species[[s_i]], '_count'),
          row[[i]],
          timestep
        )
      }
      counts <- counts + row
    }
  }
}

#' @title Step mosquito solver
#' @description calculates total_M per species and updates the vector ode
#'
#' @param solvers for each species
#' @noRd
create_solver_stepping_process <- function(solvers, parameters) {
  function(timestep) {
    for (i in seq_along(solvers)) {
      if (parameters$species_proportions[[i]] > 0) {
        solvers[[i]]$step()
      }
    }
  }
}

Solver <- R6::R6Class(
  'Solver',
  private = list(
    .solver = NULL
  ),
  public = list(
    initialize = function(solver) {
      private$.solver <- solver
    },
    step = function() {
      solver_step(private$.solver)
    },
    get_states = function() {
      solver_get_states(private$.solver)
    },

    # This is the same as `get_states`, just exposed under the interface that
    # is expected of stateful objects.
    save_state = function() {
      solver_get_states(private$.solver)
    },
    restore_state = function(t, state) {
      solver_set_states(private$.solver, t, state)
    }
  )
)

AquaticMosquitoModel <- R6::R6Class(
  'AquaticMosquitoModel',
  public = list(
    .model = NULL,
    initialize = function(model) {
      self$.model <- model
    },

    # The aquatic mosquito model doesn't have any state to save or restore (the
    # state of the ODE is stored separately). We still provide these methods to
    # conform to the expected interface.
    save_state = function() { NULL },
    restore_state = function(t, state) { }
  )
)

AdultMosquitoModel <- R6::R6Class(
  'AdultMosquitoModel',
  public = list(
    .model = NULL,
    initialize = function(model) {
      self$.model <- model
    },
    save_state = function() {
      adult_mosquito_model_save_state(self$.model)
    },
    restore_state = function(t, state) {
      adult_mosquito_model_restore_state(self$.model, state)
    }
  )
)

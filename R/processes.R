#' @title Define model processes
#' @description
#' create_processes, defines the functions which describe how each individual's
#' variables change over time.
#'
#' @param renderer a renderer object
#' @param variables a list of variables in the model
#' @param events a list of events in the model
#' @param parameters a list of model parameters
#' @param models a list of vector models, one for each species
#' @param solvers a list of ode solvers, one for each species
#' @param correlations the intervention correlations object
#' @noRd
create_processes <- function(
  renderer,
  variables,
  events,
  parameters,
  models,
  solvers,
  correlations
  ) {
  # ========
  # Immunity
  # ========
  processes <- list(
    # Maternal immunity
    create_exponential_decay_process(variables$icm, parameters$rm),
    create_exponential_decay_process(variables$ivm, parameters$rvm),
    # Blood immunity
    create_exponential_decay_process(variables$ib, parameters$rb),
    # Acquired immunity
    create_exponential_decay_process(variables$ica, parameters$rc),
    create_exponential_decay_process(variables$iva, parameters$rva),
    create_exponential_decay_process(variables$id, parameters$rid)
  )

  if (parameters$individual_mosquitoes) {
    processes <- c(
      processes,
      create_mosquito_emergence_process(
        solvers,
        variables$mosquito_state,
        variables$species,
        parameters$species,
        parameters$dpl
      )
    )
  }

  # ==============================
  # Biting and mortality processes
  # ==============================
  # schedule infections for humans and set last_boosted_*
  # move mosquitoes into incubating state
  # kill mosquitoes caught in vector control
  init_eir <- lapply(
    seq_along(parameters$species),
    function(species) {
      LaggedValue$new(
        parameters$de + 1,
        calculate_eir(
          species,
          solvers,
          variables,
          parameters,
          0
        )
      )
    }
  )

  processes <- c(
    processes,
    create_biting_process(
      renderer,
      solvers,
      models,
      variables,
      events,
      parameters,
      LaggedValue$new(
        parameters$delay_gam + 2,
        parameters$init_foim
      ), # lagged FOIM
      init_eir # lagged EIR
    ),
    create_mortality_process(variables, events, renderer, parameters),
    create_progression_process(
      events$asymptomatic_progression,
      variables$state,
      'D',
      parameters$dd
    ),
    create_progression_process(
      events$subpatent_progression,
      variables$state,
      'A',
      parameters$da
    ),
    create_progression_process(
      events$recovery,
      variables$state,
      'U',
      parameters$du
    ),
    create_progression_process(
      events$recovery,
      variables$state,
      'Tr',
      parameters$dt
    )
  )

  # ===============
  # ODE integration
  # ===============
  processes <- c(
    processes,
    create_solver_stepping_process(solvers)
  )

  # =========
  # Rendering
  # =========
  processes <- c(
    processes,
    individual::categorical_count_renderer_process(
      renderer,
      variables$state,
      c('S', 'A', 'D', 'U', 'Tr')
    ),
    create_variable_mean_renderer_process(
      renderer,
      c('ica', 'icm', 'ib', 'id'),
      variables[c('ica', 'icm', 'ib', 'id')]
    ),
    create_prevelance_renderer(
      variables$state,
      variables$birth,
      variables$is_severe,
      variables$id,
      parameters,
      renderer
    ),
    create_ode_rendering_process(renderer, solvers, parameters)
  )

  if (parameters$individual_mosquitoes) {
    processes <- c(
      processes,
      individual::categorical_count_renderer_process(
        renderer,
        variables$mosquito_state,
        c('Sm', 'Pm', 'Im')
      ),
      create_total_M_renderer_individual(
        variables$mosquito_state,
        variables$species,
        renderer,
        parameters
      )
    )
  } else {
    processes <- c(
      processes,
      create_total_M_renderer_compartmental(
        renderer,
        solvers
      )
    )
  }

  # ======================
  # Intervention processes
  # ======================

  if (parameters$bednets) {
    processes <- c(
      processes,
      distribute_nets(
        variables,
        events$throw_away_net,
        parameters,
        correlations
      )
    )
  }

  if (parameters$spraying) {
    processes <- c(
      processes,
      indoor_spraying(variables$spray_time, parameters, correlations)
    )
  }

  processes
}

# =================
# Utility functions
# =================

#' @title Exponentially decaying variables
#' @description
#' create_exponential_decay_process generates a process function
#' that reduces the value of a variable at an exponential rate
#'
#' @param variable the variable to update
#' @param rate the exponential rate
#' @noRd
create_exponential_decay_process <- function(variable, rate) {
  decay_rate <- exp(-1/rate)
  function(timestep) variable$queue_update(variable$get_values() * decay_rate)
}

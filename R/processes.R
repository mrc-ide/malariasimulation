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
#' @param lagged_eir a list of list of LaggedValue objects for EIR for each
#' population and species in the simulation
#' @param lagged_infectivity a list of LaggedValue objects for FOIM for each population
#' in the simulation
#' @param timesteps Number of timesteps
#' @param mixing a vector of mixing coefficients for the lagged transmission
#' values (default: 1)
#' @param mixing_index an index for this population's position in the
#' lagged transmission lists (default: 1)
#' @noRd
create_processes <- function(
    renderer,
    variables,
    events,
    parameters,
    models,
    solvers,
    correlations,
    lagged_eir,
    lagged_infectivity,
    timesteps,
    mixing = 1,
    mixing_index = 1
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
  processes <- c(
    processes,
    create_biting_process(
      renderer,
      solvers,
      models,
      variables,
      events,
      parameters,
      lagged_infectivity,
      lagged_eir,
      mixing,
      mixing_index
    ),
    create_asymptomatic_progression_process(
      variables$state,
      parameters$dd,
      variables,
      parameters
    ),
    create_progression_process(
      variables$state,
      'A',
      'U',
      parameters$da,
      variables$infectivity,
      parameters$cu
    ),
    create_progression_process(
      variables$state,
      'U',
      'S',
      parameters$du,
      variables$infectivity,
      0
    )
  )

  # =======================
  # Antimalarial Resistance
  # =======================
  # Add an a new process which governs the transition from Tr to S when
  # antimalarial resistance is simulated. The rate of transition switches
  # from a parameter to a variable when antimalarial resistance == TRUE.

  # Assign the dt input to a separate object with the default single parameter value:
  dt_input <- parameters$dt

  # If antimalarial resistance is switched on, assign dt variable values to the
  if(parameters$antimalarial_resistance) {
    dt_input <- variables$dt
  }

  # Create the progression process for Tr --> S specifying dt_input as the rate:
  processes <- c(
    processes,
    create_progression_process(
      variables$state,
      'Tr',
      'S',
      dt_input,
      variables$infectivity,
      0
    )
  )

  # ===============
  # ODE integration
  # ===============
  processes <- c(
    processes,
    create_solver_stepping_process(solvers, parameters)
  )

  # =========
  # RTS,S EPI
  # =========
  if (!is.null(parameters$pev_epi_coverage)) {
    processes <- c(
      processes,
      create_epi_pev_process(
        variables,
        events,
        parameters,
        correlations,
        parameters$pev_epi_coverage,
        parameters$pev_epi_timesteps
      )
    )
  }

  # =========
  # PMC
  # =========
  if(!is.null(parameters$pmc_coverages)){
    processes <- c(
      processes,
      create_pmc_process(
        variables,
        events,
        parameters,
        renderer,
        correlations,
        parameters$pmc_coverages,
        parameters$pmc_timesteps,
        parameters$pmc_drug
      )
    )
  }

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
      c('ica', 'icm', 'ib', 'id', 'iva', 'ivm'),
      variables[c('ica', 'icm', 'ib', 'id', 'iva', 'ivm')]
    ),
    create_prevelance_renderer(
      variables$state,
      variables$birth,
      variables$id,
      parameters,
      renderer
    ),
    create_age_group_renderer(
      variables$birth,
      parameters,
      renderer
    ),
    create_compartmental_rendering_process(renderer, solvers, parameters)
  )

  if (parameters$individual_mosquitoes) {
    processes <- c(
      processes,
      create_vector_count_renderer_individual(
        variables$mosquito_state,
        variables$species,
        variables$mosquito_state,
        renderer,
        parameters
      )
    )
  } else {
    processes <- c(
      processes,
      create_total_M_renderer_compartmental(
        renderer,
        solvers,
        parameters
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
      ),
      net_usage_renderer(variables$net_time, renderer)
    )
  }

  if (parameters$spraying) {
    processes <- c(
      processes,
      indoor_spraying(variables$spray_time, parameters, correlations),
      spray_renderer(variables$spray_time, renderer)
    )
  }

  # ======================
  # Progress bar process
  # ======================
  if (parameters$progress_bar){
    processes <- c(
      processes,
      create_progress_process(timesteps)
    )
  }

  # Mortality step
  processes <- c(
    processes,
    create_mortality_process(variables, events, renderer, parameters)
  )

  # Combined interventions
  # pev - bednets
  processes <- c(
    processes,
    create_combined_intervention_rendering_process(
      'pev',
      variables$last_pev_timestep,
      'bednets',
      variables$net_time,
      365,
      renderer
    )
  )

  # smc - bednets
  processes <- c(
    processes,
    create_combined_intervention_rendering_process(
      'smc',
      variables$drug_time,
      'bednets',
      variables$net_time,
      365,
      renderer
    )
  )

  # pev - smc
  processes <- c(
    processes,
    create_combined_intervention_rendering_process(
      'pev',
      variables$last_pev_timestep,
      'smc',
      variables$drug_time,
      365,
      renderer
    )
  )

  # irs - bednet
  processes <- c(
    processes,
    create_combined_intervention_rendering_process(
      'spraying',
      variables$spray_time,
      'bednets',
      variables$net_time,
      365,
      renderer
    )
  )

  # irs - smc
  processes <- c(
    processes,
    create_combined_intervention_rendering_process(
      'spraying',
      variables$spray_time,
      'smc',
      variables$drug_time,
      365,
      renderer
    )
  )

  # irs - pev
  processes <- c(
    processes,
    create_combined_intervention_rendering_process(
      'spraying',
      variables$spray_time,
      'pev',
      variables$last_pev_timestep,
      365,
      renderer
    )
  )

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
  stopifnot(inherits(variable, "DoubleVariable"))
  decay_rate <- exp(-1/rate)
  exponential_process_cpp(variable$.variable, decay_rate)
}

#' @title Create and initialise lagged_infectivity object
#'
#' @param variables model variables for initialisation
#' @param parameters model parameters
#' @noRd
create_lagged_infectivity <- function(variables, parameters) {
  age <- get_age(variables$birth$get_values(), 0)
  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(psi, variables$zeta$get_values())
  init_infectivity <- sum(.pi * variables$infectivity$get_values())
  LaggedValue$new(
    max_lag = parameters$delay_gam + 2,
    default = init_infectivity
  )
}

#' @title Create and initialise a list of lagged_eir objects per species
#'
#' @param variables model variables for initialisation
#' @param solvers model differential equation solvers
#' @param parameters model parameters
#' @noRd
create_lagged_eir <- function(variables, solvers, parameters) {
  lapply(
    seq_along(parameters$species),
    function(species) {
      LaggedValue$new(
        max_lag = parameters$de + 1,
        default = calculate_eir(
          species,
          solvers,
          variables,
          parameters,
          0
        )
      )
    }
  )
}

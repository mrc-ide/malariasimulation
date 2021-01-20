#' @title Define model processes
#' @description
#' create_processes, defines the functions which describe how each individual's
#' states and variables change over time.
#'
#' It lists processes from `infection.R`, `mosquito_emergence.R` and
#' `mortality.R`; and then exposes them to the model
#' @param variables a list of variables in the model
#' @param events a list of events in the model
#' @param parameters a list of model parameters
#' @param odes a list of vector ode models for each species
create_processes <- function(
  variables,
  events,
  parameters,
  odes
  ) {
  processes <- list(
    # ========
    # Immunity
    # ========

    # Maternal immunity
    create_exponential_decay_process(variables$icm, parameters$rm),
    create_exponential_decay_process(variables$ivm, parameters$rvm),
    # Blood immunity
    create_exponential_decay_process(variables$ib, parameters$rb),
    # Acquired immunity
    create_exponential_decay_process(variables$ica, parameters$rc),
    create_exponential_decay_process(variables$iva, parameters$rva),
    create_exponential_decay_process(variables$id, parameters$rid),

    create_mosquito_emergence_process_cpp(
      odes,
      states$Unborn$name,
      states$Sm$name,
      variables$mosquito_species$.variable,
      parameters$dpl
    ),

    # ==============
    # Biting process
    # ==============
    # schedule infections for humans and set last_boosted_*
    # move mosquitoes into incubating state
    # kill mosquitoes caught in vector control
    create_biting_process(variables, events, parameters),

    create_mortality_process(),

    # ===============
    # ODE integration
    # ===============
    create_ode_stepping_process_cpp(
      odes,
      variables$mosquito_state,
      c('Sm', 'Pm', 'Im'),
      variables$mosquito_species
    ),

    # Rendering processes
    individual::state_count_renderer_process(
      variables$state,
      c('S', 'A', 'D', 'U', 'Tr')
    ),
    individual::variable_mean_renderer_process(
      individuals$human$name,
      c(
        variables$ica$name,
        variables$icm$name,
        variables$ib$name,
        variables$id$name
      )
    ),
    create_prevelance_renderer(
      individuals$human,
      variables$state,
      variables$birth,
      variables$is_severe
    ),
    individual::state_count_renderer_process(
      variables$mosquito_state,
      c('Sm', 'Pm', 'Im')
    ),
    
    create_ode_rendering_process(odes)
  )

  if (parameters$bednets) {
    processes <- c(
      processes,
      distribute_nets(
        variables,
        events$throw_away_net,
        parameters
      )
    )
  }

  if (parameters$spraying) {
    processes <- c(
      processes,
      indoor_spraying(variables$spray_time, parameters)
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
create_exponential_decay_process <- function(variable, rate) {
  decay_rate <- exp(-1/rate)
  function(timestep) variable$queue_update(variable$get_values() * decay_rate)
}

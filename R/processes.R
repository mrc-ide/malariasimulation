#' @title Define model processes
#' @description
#' create_processes, defines the functions which describe how each individual's
#' variables change over time.
#'
#' @param renderer a renderer object
#' @param variables a list of variables in the model
#' @param events a list of events in the model
#' @param parameters a list of model parameters
#' @param odes a list of vector ode models for each species
create_processes <- function(
  renderer,
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
      variables$mosquito_state$.variable,
      variables$species$.variable,
      parameters$species,
      parameters$dpl
    ),

    # ==============
    # Biting process
    # ==============
    # schedule infections for humans and set last_boosted_*
    # move mosquitoes into incubating state
    # kill mosquitoes caught in vector control
    create_biting_process(renderer, variables, events, parameters),

    create_mortality_process(variables, events, renderer, parameters)

    ## ===============
    ## ODE integration
    ## ===============
    #create_ode_stepping_process_cpp(
      #odes,
      #variables$mosquito_state$.variable,
      #variables$species$.variable,
      #parameters$species
    #),

    ## Rendering processes
    #individual::categorical_count_renderer_process(
      #renderer,
      #variables$state,
      #c('S', 'A', 'D', 'U', 'Tr')
    #),
    #create_variable_mean_renderer_process(
      #renderer,
      #c('ica', 'icm', 'ib', 'id'),
      #variables[c('ica', 'icm', 'ib', 'id')]
    #),
    #create_prevelance_renderer(
      #variables$state,
      #variables$birth,
      #variables$is_severe,
      #parameters,
      #renderer
    #),
    #individual::categorical_count_renderer_process(
      #renderer,
      #variables$mosquito_state,
      #c('Sm', 'Pm', 'Im')
    #),
    
    #create_ode_rendering_process(renderer, odes)
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

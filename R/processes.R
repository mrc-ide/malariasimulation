#' @title Define model processes
#' @description
#' create_processes, defines the functions which describe how each individual's
#' states and variables change over time.
#'
#' It lists processes from `infection.R`, `mosquito_emergence.R` and
#' `mortality.R`; and then exposes them to the model
#' @param individuals a list of individuals in the model
#' @param states a list of states in the model
#' @param variables a list of variables in the model
#' @param events a list of events in the model
#' @param parameters a list of model parameters
#' @param odes a list of vector models (if vector_ode == TRUE)
create_processes <- function(
  individuals,
  states,
  variables,
  events,
  parameters,
  odes=NULL
  ) {
  processes <- list(
    # ========
    # Immunity
    # ========

    # Maternal immunity
    create_exponential_decay_process(individuals$human, variables$icm, parameters$rm),
    create_exponential_decay_process(individuals$human, variables$ivm, parameters$rvm),
    # Blood immunity
    create_exponential_decay_process(individuals$human, variables$ib, parameters$rb),
    # Acquired immunity
    create_exponential_decay_process(individuals$human, variables$ica, parameters$rc),
    create_exponential_decay_process(individuals$human, variables$iva, parameters$rva),
    create_exponential_decay_process(individuals$human, variables$id, parameters$rid),

    # =================
    # State transitions
    # =================
    create_asymptomatic_progression_process(
      individuals$human,
      states,
      variables,
      1. - exp(-1./parameters$dd)
    ),
    create_progression_process(
      individuals$human,
      states$A,
      states$U,
      1. - exp(-1./parameters$da),
      variables$infectivity,
      parameters$cu
    ),
    create_progression_process(
      individuals$human,
      states$U,
      states$S,
      1. - exp(-1./parameters$du),
      variables$infectivity,
      0
    ),
    create_progression_process(
      individuals$human,
      states$Tr,
      states$S,
      1. - exp(-1./parameters$dt),
      variables$infectivity,
      0
    ),

    # schedule infections for humans and set last_boosted_*
    create_infection_process(
      individuals,
      states,
      variables,
      events,
      odes
    ),

    create_mortality_process(
      individuals$human,
      states$D,
      states$Tr,
      variables,
      events
    ),

    # Rendering processes
    individual::state_count_renderer_process(
      individuals$human$name,
      c(
        states$S$name,
        states$A$name,
        states$D$name,
        states$U$name,
        states$Tr$name
      )
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
      states$D,
      states$A,
      variables$birth
    )
  )

  if (!parameters$vector_ode) {
    processes <- c(
      processes,
      # ==================
      # Mosquito Processes
      # ==================
      create_mosquito_infection_process(
        individuals$mosquito,
        individuals$human,
        states,
        variables,
        events$mosquito_infection
      ),

      # Eggs laid
      create_egg_laying_process_cpp(
        individuals$mosquito$name,
        states$Sm$name,
        states$Pm$name,
        states$Im$name,
        states$Unborn$name,
        states$E$name
      ),

      individual::fixed_probability_state_change_process(
        individuals$mosquito$name,
        states$E$name,
        states$L$name,
        1. - exp(-1./parameters$del)
      ),

      # Death of larvae
      create_larval_death_process_cpp(
        individuals$mosquito$name,
        states$E$name,
        states$L$name,
        states$Unborn$name,
        calculate_carrying_capacity(parameters),
        calculate_R_bar(parameters)
      ),

      individual::fixed_probability_state_change_process(
        individuals$mosquito$name,
        states$L$name,
        states$P$name,
        1. - exp(-1./parameters$dl)
      ),

      individual::fixed_probability_state_change_process(
        individuals$mosquito$name,
        states$P$name,
        states$Sm$name,
        .5 * (1. - exp(-1./parameters$dpl))
      ),

      # Death of pupals
      individual::fixed_probability_state_change_process(
        individuals$mosquito$name,
        states$P$name,
        states$Unborn$name,
        parameters$mup
      ),

      # Natural death of females
      individual::fixed_probability_state_change_process(
        individuals$mosquito$name,
        states$Sm$name,
        states$Unborn$name,
        parameters$mum
      ),
      individual::fixed_probability_state_change_process(
        individuals$mosquito$name,
        states$Pm$name,
        states$Unborn$name,
        parameters$mum
      ),
      individual::fixed_probability_state_change_process(
        individuals$mosquito$name,
        states$Im$name,
        states$Unborn$name,
        parameters$mum
      ),

      # Rendering processes
      individual::state_count_renderer_process(
        individuals$mosquito$name,
        c(
          states$E$name,
          states$L$name,
          states$P$name,
          states$Sm$name,
          states$Pm$name,
          states$Im$name
        )
      )
    )
  } else {
    processes <- c(
      processes,
      create_ode_stepping_process(odes, individuals$human, states, variables),
      create_ode_rendering_process(odes)
    )
  }
}

#' @title Define event based processes
#' @description defines processes for events that can be scheduled in the future
#'
#' @param individuals a list of individuals in the model
#' @param states a list of states in the model
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
create_event_based_processes <- function(individuals, states, variables, events, parameters) {
  events$infection$add_listener(
    individual::update_state_listener(individuals$human$name, states$D$name)
  )
  events$infection$add_listener(
    function(api, target) {
      api$queue_variable_update(
        individuals$human,
        variables$infectivity,
        parameters$cd,
        target
      )
    }
  )
  events$asymptomatic_infection$add_listener(
    individual::update_state_listener(individuals$human$name, states$A$name)
  )
  events$asymptomatic_infection$add_listener(
    function(api, target) {
      new_infectivity <- asymptomatic_infectivity(
        get_age(
          api$get_variable(individuals$human, variables$birth, target),
          api$get_timestep()
        ),
        api$get_variable(individuals$human, variables$id, target),
        api$get_parameters()
      )
      api$queue_variable_update(
        individuals$human,
        variables$infectivity,
        new_infectivity,
        target
      )
    }
  )

  if (!parameters$vector_ode) {
    events$mosquito_infection$add_listener(
      individual::update_state_listener(individuals$mosquito$name, states$Im$name)
    )
  }
}

# =================
# Utility functions
# =================

#' @title Exponentially decaying variables
#' @description
#' create_exponential_decay_process generates a process function
#' that reduces the value of a variable at an exponential rate
#'
#' @param individual an individual
#' @param variable the variable to update
#' @param rate the exponential rate
create_exponential_decay_process <- function(individual, variable, rate) {
  decay_rate <- exp(-1/rate)
  function(api) {
    i <- api$get_variable(individual, variable)
    api$queue_variable_update(individual, variable, i * decay_rate)
  }
}

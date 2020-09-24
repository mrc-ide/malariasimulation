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
#' @param odes a list of vector ode models for each species
#' @param mda_events a list of lists of events for each MDA
create_processes <- function(
  individuals,
  states,
  variables,
  events,
  parameters,
  odes,
  mda_events
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
      events
    ),

    create_mortality_process_cpp(),

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
      variables$birth,
      variables$is_severe
    ),

    # ==================
    # Mosquito Processes
    # ==================
    create_mosquito_infection_process_cpp(
      individuals$mosquito$name,
      individuals$human$name,
      c(states$Sm$name, states$Pm$name),
      c(
        variables$birth$name,
        variables$zeta$name,
        variables$infectivity$name,
        variables$mosquito_variety$name
      )
    ),

    create_mosquito_emergence_process_cpp(
      individuals$mosquito$name,
      odes,
      states$Unborn$name,
      states$Sm$name,
      variables$mosquito_variety$name,
      parameters$dpl
    ),

    # Infection after incubation
    individual::fixed_probability_state_change_process(
      individuals$mosquito$name,
      states$Pm$name,
      states$Im$name,
      1. - exp(-1./parameters$dem)
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
        states$Sm$name,
        states$Pm$name,
        states$Im$name
      )
    ),
    create_ode_stepping_process_cpp(
      odes,
      individuals$mosquito$name,
      c(states$Sm$name, states$Pm$name, states$Im$name),
      variables$mosquito_variety$name
    ),
    create_ode_rendering_process(odes)
  )
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
  events$clinical_infection$add_listener(
    individual::update_state_listener(individuals$human$name, states$D$name)
  )
  events$clinical_infection$add_listener(
    function(api, target) {
      if (length(target) > 0) {
        api$queue_variable_update(
          individuals$human,
          variables$infectivity,
          parameters$cd,
          target
        )
      }
    }
  )
  events$asymptomatic_infection$add_listener(
    individual::update_state_listener(individuals$human$name, states$A$name)
  )
  events$asymptomatic_infection$add_listener(
    function(api, target) {
      if (length(target) > 0) {
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
    }
  )

  events$infection$add_listener(
    create_incidence_renderer(
      individuals$human,
      variables$birth,
      variables$is_severe
    )
  )

  if (parameters$rtss == 1) {
    events$rtss_vaccination$add_listener(
      create_rtss_vaccination_listener(individuals$human, variables, events, parameters)
    )
    events$rtss_booster$add_listener(
      create_rtss_booster_listener(individuals$human, variables, events, parameters)
    )
  }

  if (parameters$mda == 1) {
    mda_listeners <- create_mda_listeners(
      individuals$human,
      states,
      variables,
      events$mda_administer,
      parameters$mda_drug,
      parameters$mda_end,
      parameters$mda_frequency,
      parameters$mda_min_age,
      parameters$mda_max_age,
      parameters$mda_coverage
    )
    events$mda_enrollment$add_listener(mda_listeners$enrollment_listener)
    events$mda_administer$add_listener(mda_listeners$administer_listener)
  }

  if (parameters$smc == 1) {
    smc_listeners <- create_mda_listeners(
      individuals$human,
      states,
      variables,
      events$smc_administer,
      parameters$smc_drug,
      parameters$smc_end,
      parameters$smc_frequency,
      parameters$smc_min_age,
      parameters$smc_max_age,
      parameters$smc_coverage
    )
    events$smc_enrollment$add_listener(smc_listeners$enrollment_listener)
    events$smc_administer$add_listener(smc_listeners$administer_listener)
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

create_setup_process <- function(events) {
  function(api) {
    parameters <- api$get_parameters()
    if (parameters$rtss) {
      api$schedule(events$rtss_vaccination, c(1), parameters$rtss_start)
    }
    if (parameters$mda) {
      api$schedule(events$mda_enrollment, c(1), parameters$mda_start)
    }
    if (parameters$smc) {
      api$schedule(events$smc_enrollment, c(1), parameters$smc_start)
    }
  }
}

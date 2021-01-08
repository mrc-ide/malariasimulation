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
create_processes <- function(
  individuals,
  states,
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
    create_exponential_decay_process(individuals$human, variables$icm, parameters$rm),
    create_exponential_decay_process(individuals$human, variables$ivm, parameters$rvm),
    # Blood immunity
    create_exponential_decay_process(individuals$human, variables$ib, parameters$rb),
    # Acquired immunity
    create_exponential_decay_process(individuals$human, variables$ica, parameters$rc),
    create_exponential_decay_process(individuals$human, variables$iva, parameters$rva),
    create_exponential_decay_process(individuals$human, variables$id, parameters$rid),

    create_mosquito_emergence_process_cpp(
      individuals$mosquito$name,
      odes,
      states$Unborn$name,
      states$Sm$name,
      variables$mosquito_variety$name,
      parameters$dpl
    ),

    # ==============
    # Biting process
    # ==============
    # schedule infections for humans and set last_boosted_*
    # move mosquitoes into incubating state
    # kill mosquitoes caught in vector control
    create_biting_process(
      individuals,
      states,
      variables,
      events
    ),

    create_mortality_process_cpp(),

    # ===============
    # ODE integration
    # ===============
    create_ode_stepping_process_cpp(
      odes,
      individuals$mosquito$name,
      c(states$Sm$name, states$Pm$name, states$Im$name),
      variables$mosquito_variety$name
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
      variables$birth,
      variables$is_severe
    ),
    individual::state_count_renderer_process(
      individuals$mosquito$name,
      c(
        states$Sm$name,
        states$Pm$name,
        states$Im$name
      )
    ),
    
    create_ode_rendering_process(odes)
  )

  if (parameters$bednets) {
    processes <- c(
      processes,
      distribute_nets(
        individuals$human,
        variables,
        events$throw_away_net,
        parameters
      )
    )
  }

  if (parameters$spraying) {
    processes <- c(
      processes,
      indoor_spraying(individuals$human, variables$spray_time, parameters)
    )
  }

  processes
}

#' @title Define event based processes
#' @description defines processes for events that can be scheduled in the future
#'
#' @param individuals a list of individuals in the model
#' @param states a list of states in the model
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
create_event_based_processes <- function(
  individuals,
  states,
  variables,
  events,
  parameters,
  correlations
  ) {

  # =============
  # State updates
  # =============
  # When infection events fire, update the corresponding states and infectivity
  # variables

  # Infection events
  events$clinical_infection$add_listener(
    create_infection_update_listener(
      individuals$human,
      states$D,
      variables$infectivity,
      parameters$cd
    )
  )
  events$asymptomatic_infection$add_listener(
    create_asymptomatic_update_listener(
      individuals$human,
      states,
      variables
    )
  )

  # Recovery events
  events$subpatent_infection$add_listener(
    create_infection_update_listener(
      individuals$human,
      states$U,
      variables$infectivity,
      parameters$cu
    )
  )
  events$recovery$add_listener(
    create_infection_update_listener(
      individuals$human,
      states$S,
      variables$infectivity,
      0
    )
  )

  # ===========
  # Progression
  # ===========
  # When infection events fire, schedule the next stages of infection

  events$clinical_infection$add_listener(
    create_progression_listener(
      events$asymptomatic_infection,
      parameters$dd
    )
  )
  events$asymptomatic_infection$add_listener(
    create_progression_listener(
      events$subpatent_infection,
      parameters$da
    )
  )
  events$subpatent_infection$add_listener(
    create_progression_listener(
      events$recovery,
      parameters$du
    )
  )

  events$infection$add_listener(
    create_incidence_renderer(
      individuals$human,
      variables$birth,
      variables$is_severe
    )
  )

  events$mosquito_infection$add_listener(
    individual::update_state_listener(
      individuals$mosquito$name,
      states$Im$name
    )
  )

  if (parameters$bednets == 1) {
    events$throw_away_net$add_listener(
      throw_away_nets(individuals$human, variables)
    )
  }

  if (parameters$rtss == 1) {
    events$rtss_vaccination$add_listener(
      create_rtss_vaccination_listener(
        individuals$human,
        variables,
        events,
        parameters,
        correlations
      )
    )
    events$rtss_booster$add_listener(
      create_rtss_booster_listener(
        individuals$human,
        variables,
        events,
        parameters
      )
    )
  }

  if (parameters$mda == 1) {
    events$mda_administer$add_listener(create_mda_listeners(
      individuals$human,
      states,
      variables,
      events$mda_administer,
      parameters$mda_drug,
      parameters$mda_end,
      parameters$mda_frequency,
      parameters$mda_min_age,
      parameters$mda_max_age,
      parameters$mda_coverage,
      correlations,
      'mda'
    ))
  }

  if (parameters$smc == 1) {
    events$smc_administer$add_listener(create_mda_listeners(
      individuals$human,
      states,
      variables,
      events$smc_administer,
      parameters$smc_drug,
      parameters$smc_end,
      parameters$smc_frequency,
      parameters$smc_min_age,
      parameters$smc_max_age,
      parameters$smc_coverage,
      correlations,
      'smc'
    ))
  }

  events$tbv_vaccination$add_listener(
    function(api, target) {
      timestep <- api$get_timestep()
      target <- which(trunc(get_age(
        api$get_variable(individuals$human, variables$birth),
        timestep
      ) / 365) %in% parameters$tbv_ages)
      to_vaccinate <- which(sample_intervention(
        target,
        'tbv',
        parameters$tbv_coverage,
        correlations
      ))
      api$render('n_vaccinated_tbv', length(target))
      if (length(to_vaccinate) > 0) {
        api$queue_variable_update(
          individuals$human,
          variables$tbv_vaccinated,
          timestep,
          to_vaccinate
        )
      }
      if (timestep + parameters$tbv_frequency <= parameters$tbv_end) {
        api$schedule(events$tbv_vaccination, c(1), parameters$tbv_frequency)
      }
    }
  )
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

create_setup_process <- function(individuals, states, events) {
  function(api) {
    parameters <- api$get_parameters()
    # Initialise malaria progression
    initialise_progression(
      api,
      events$asymptomatic_infection,
      individuals$human,
      states$D,
      parameters$de
    )
    initialise_progression(
      api,
      events$subpatent_infection,
      individuals$human,
      states$A,
      parameters$da
    )
    initialise_progression(
      api,
      events$recovery,
      individuals$human,
      states$U,
      parameters$du
    )
    initialise_progression(
      api,
      events$recovery,
      individuals$human,
      states$Tr,
      parameters$dt
    )
    api$schedule(
      events$mosquito_infection,
      api$get_state(individuals$mosquito, states$Pm),
      parameters$dem
    )

    # Initialise interventions
    if (parameters$rtss) {
      api$schedule(events$rtss_vaccination, c(1), parameters$rtss_start)
    }
    if (parameters$mda) {
      api$schedule(events$mda_administer, c(1), parameters$mda_start)
    }
    if (parameters$smc) {
      api$schedule(events$smc_administer, c(1), parameters$smc_start)
    }
    if (parameters$tbv) {
      api$schedule(events$tbv_vaccination, c(1), parameters$tbv_start)
    }
  }
}

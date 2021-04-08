create_events <- function(parameters) {
  events <- list(
    # Human infection events
    clinical_infection = individual::TargetedEvent$new(parameters$human_population),
    asymptomatic_infection = individual::TargetedEvent$new(parameters$human_population),
    # either clinical or asym infection
    infection = individual::TargetedEvent$new(parameters$human_population), 
    subpatent_infection = individual::TargetedEvent$new(parameters$human_population),
    recovery = individual::TargetedEvent$new(parameters$human_population),

    # Vaccination events
    rtss_vaccination = individual::Event$new(),
    rtss_booster = individual::TargetedEvent$new(parameters$human_population),
    tbv_vaccination = individual::Event$new(),

    # MDA events
    mda_administer = individual::Event$new(),
    smc_administer = individual::Event$new(),

    # Bednet events
    throw_away_net = individual::TargetedEvent$new(parameters$human_population)
  )
  if (parameters$individual_mosquitoes) {
    events <- c(
      events,
      # Mosquito events
      mosquito_infection = individual::TargetedEvent$new(parameters$mosquito_limit),
      mosquito_death = individual::TargetedEvent$new(parameters$mosquito_limit)
    )
  }

  events
}

initialise_events <- function(events, variables, parameters) {
  initialise_progression(
    events$asymptomatic_infection,
    variables$state,
    'D',
    parameters$de
  )
  initialise_progression(
    events$subpatent_infection,
    variables$state,
    'A',
    parameters$da
  )
  initialise_progression(
    events$recovery,
    variables$state,
    'U',
    parameters$du
  )
  initialise_progression(
    events$recovery,
    variables$state,
    'Tr',
    parameters$dt
  )

  if (parameters$individual_mosquitoes) {
    incubating <- variables$mosquito_state$get_index_of('Pm')
    events$mosquito_infection$schedule(
      incubating,
      log_uniform(incubating$size(), parameters$dem)
    )
  }

  # Initialise interventions
  if (parameters$rtss) {
    events$rtss_vaccination$schedule(parameters$rtss_timesteps[[1]] - 1)
  }
  if (parameters$mda) {
    events$mda_administer$schedule(parameters$mda_timesteps[[1]] - 1)
  }
  if (parameters$smc) {
    events$smc_administer$schedule(parameters$smc_timesteps[[1]] - 1)
  }
  if (parameters$tbv) {
    events$tbv_vaccination$schedule(parameters$tbv_timesteps[[1]] - 1)
  }
}

#' @title Define event based processes
#' @description defines processes for events that can be scheduled in the future
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @param renderer the model rendering object
#' @noRd
attach_event_listeners <- function(
  events,
  variables,
  parameters,
  correlations,
  renderer
  ) {

  # =============
  # State updates
  # =============
  # When infection events fire, update the corresponding states and infectivity
  # variables

  # Infection events
  events$clinical_infection$add_listener(
    create_infection_update_listener(
      variables$state,
      'D',
      variables$infectivity,
      parameters$cd
    )
  )

  events$asymptomatic_infection$add_listener(
    create_asymptomatic_update_listener(
      variables,
      parameters
    )
  )

  # Recovery events
  events$subpatent_infection$add_listener(
    create_infection_update_listener(
      variables$state,
      'U',
      variables$infectivity,
      parameters$cu
    )
  )
  events$recovery$add_listener(
    create_infection_update_listener(
      variables$state,
      'S',
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
      variables$birth,
      variables$is_severe,
      parameters,
      renderer
    )
  )

  if (parameters$individual_mosquitoes) {
    events$mosquito_infection$add_listener(
      individual::update_category_listener(variables$mosquito_state, 'Im')
    )

    events$mosquito_death$add_listener(
      function(timestep, target) {
        variables$mosquito_state$queue_update('NonExistent', target)
        events$mosquito_infection$clear_schedule(target)
        renderer$render('mosquito_deaths', target$size(), timestep)
      }
    )
  }

  if (parameters$bednets == 1) {
    events$throw_away_net$add_listener(
      throw_away_nets(variables)
    )
  }

  if (parameters$rtss == 1) {
    events$rtss_vaccination$add_listener(
      create_rtss_vaccination_listener(
        variables,
        events,
        parameters,
        correlations,
        renderer
      )
    )
    events$rtss_booster$add_listener(
      create_rtss_booster_listener(
        variables,
        events,
        parameters
      )
    )
  }

  if (parameters$mda == 1) {
    events$mda_administer$add_listener(create_mda_listeners(
      variables,
      events$mda_administer,
      parameters$mda_drug,
      parameters$mda_timesteps,
      parameters$mda_coverages,
      parameters$mda_min_age,
      parameters$mda_max_age,
      correlations,
      'mda',
      parameters,
      renderer
    ))
  }

  if (parameters$smc == 1) {
    events$smc_administer$add_listener(create_mda_listeners(
      variables,
      events$smc_administer,
      parameters$smc_drug,
      parameters$smc_timesteps,
      parameters$smc_coverages,
      parameters$smc_min_age,
      parameters$smc_max_age,
      correlations,
      'smc',
      parameters,
      renderer
    ))
  }

  if (parameters$tbv == 1) {
    events$tbv_vaccination$add_listener(
      function(timestep) {
        time_index = which(parameters$tbv_timesteps == timestep)
        target <- which(trunc(get_age(
          variables$birth$get_values(),
          timestep
        ) / 365) %in% parameters$tbv_ages)
        to_vaccinate <- which(sample_intervention(
          target,
          'tbv',
          parameters$tbv_coverages[[time_index]],
          correlations
        ))
        renderer$render('n_vaccinated_tbv', length(to_vaccinate), timestep)
        if (length(to_vaccinate) > 0) {
          variables$tbv_vaccinated$queue_update(
            timestep,
            to_vaccinate
          )
        }
        if (time_index < length(parameters$tbv_timesteps)) {
          events$tbv_vaccination$schedule(
            parameters$tbv_timesteps[[time_index + 1]] - timestep
          )
        }
      }
    )
  }
}

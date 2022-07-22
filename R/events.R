create_events <- function(parameters) {
  events <- list(
    # Disease progression events
    asymptomatic_progression = individual::TargetedEvent$new(parameters$human_population),
    subpatent_progression = individual::TargetedEvent$new(parameters$human_population),
    recovery = individual::TargetedEvent$new(parameters$human_population),

    # Human infection events
    clinical_infection = individual::TargetedEvent$new(parameters$human_population),
    asymptomatic_infection = individual::TargetedEvent$new(parameters$human_population),

    # MDA events
    mda_administer = individual::Event$new(),
    smc_administer = individual::Event$new(),

    # TBV event
    tbv_vaccination = individual::Event$new(),

    # Bednet events
    throw_away_net = individual::TargetedEvent$new(parameters$human_population)
  )

  # Mass vaccination events
  if (!is.null(parameters$rtss_mass_timesteps)) {
    rtss_mass_doses <- lapply(
      seq_along(parameters$rtss_doses),
      function(.) individual::TargetedEvent$new(parameters$human_population)
    )
    rtss_mass_boosters <- lapply(
      seq_along(parameters$rtss_mass_boosters),
      function(.) individual::TargetedEvent$new(parameters$human_population)
    )
    events$rtss_mass_vaccination = individual::Event$new()
    events$rtss_mass_doses <- rtss_mass_doses
    events$rtss_mass_boosters <- rtss_mass_boosters
  }

  # EPI vaccination events
  if (!is.null(parameters$rtss_epi_coverage)) {
    rtss_epi_doses <- lapply(
      seq_along(parameters$rtss_doses),
      function(.) individual::TargetedEvent$new(parameters$human_population)
    )
    rtss_epi_boosters <- lapply(
      seq_along(parameters$rtss_epi_boosters),
      function(.) individual::TargetedEvent$new(parameters$human_population)
    )
    events$rtss_epi_doses <- rtss_epi_doses
    events$rtss_epi_boosters <- rtss_epi_boosters
  }

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
  if (parameters$individual_mosquitoes) {
    incubating <- variables$mosquito_state$get_index_of('Pm')
    events$mosquito_infection$schedule(
      incubating,
      log_uniform(incubating$size(), parameters$dem)
    )
  }

  # Initialise scheduled interventions
  if (!is.null(parameters$rtss_mass_timesteps)) {
    events$rtss_mass_vaccination$schedule(parameters$rtss_mass_timesteps[[1]] - 1)
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

  events$asymptomatic_progression$add_listener(
    create_asymptomatic_update_listener(
      variables,
      parameters
    )
  )
  events$asymptomatic_progression$add_listener(
    create_rate_listener('D', 'A', renderer)
  )
  events$asymptomatic_infection$add_listener(
    create_asymptomatic_update_listener(
      variables,
      parameters
    )
  )

  # Recovery events
  events$subpatent_progression$add_listener(
    create_infection_update_listener(
      variables$state,
      'U',
      variables$infectivity,
      parameters$cu
    )
  )

  events$subpatent_progression$add_listener(
    create_rate_listener('A', 'U', renderer)
  )

  events$recovery$add_listener(
    create_infection_update_listener(
      variables$state,
      'S',
      variables$infectivity,
      0
    )
  )
  events$recovery$add_listener(
    create_rate_listener('U', 'S', renderer)
  )

  # ===========
  # Progression
  # ===========
  # When infection events fire, schedule the next stages of infection
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

  # RTS,S listeners
  if (!is.null(events$rtss_mass_doses)) {
    # set up distribution
    events$rtss_mass_vaccination$add_listener(
      create_rtss_mass_listener(
        variables,
        events,
        parameters,
        correlations
      )
    )
    attach_rtss_dose_listeners(
      variables,
      parameters,
      events$rtss_mass_doses,
      events$rtss_mass_boosters,
      parameters$rtss_mass_boosters,
      parameters$rtss_mass_booster_coverage,
      'mass',
      renderer
    )
  }

  if (!is.null(events$rtss_epi_doses)) {
    attach_rtss_dose_listeners(
      variables,
      parameters,
      events$rtss_epi_doses,
      events$rtss_epi_boosters,
      parameters$rtss_epi_boosters,
      parameters$rtss_epi_booster_coverage,
      'epi',
      renderer
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
      create_tbv_listener(variables, events, parameters, correlations, renderer)
    )
  }
}

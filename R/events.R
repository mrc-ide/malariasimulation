create_events <- function(parameters) {
  events <- list(
    # MDA events
    mda_administer = individual::Event$new(restore=FALSE),
    smc_administer = individual::Event$new(restore=FALSE),

    # TBV event
    tbv_vaccination = individual::Event$new(restore=FALSE),

    # Bednet events
    throw_away_net = individual::TargetedEvent$new(parameters$human_population)
  )

  # Mass vaccination events
  if (!is.null(parameters$mass_pev_timesteps)) {
    mass_pev_doses <- lapply(
      seq_along(parameters$pev_doses),
      function(.) individual::TargetedEvent$new(parameters$human_population)
    )
    mass_pev_boosters <- lapply(
      seq_along(parameters$mass_pev_booster_spacing),
      function(.) individual::TargetedEvent$new(parameters$human_population)
    )
    events$mass_pev <- individual::Event$new(restore=FALSE)
    events$mass_pev_doses <- mass_pev_doses
    events$mass_pev_boosters <- mass_pev_boosters
  }

  # EPI vaccination events
  if (!is.null(parameters$pev_epi_coverage)) {
    pev_epi_doses <- lapply(
      seq_along(parameters$pev_doses),
      function(.) individual::TargetedEvent$new(parameters$human_population)
    )
    pev_epi_boosters <- lapply(
      seq_along(parameters$pev_epi_booster_spacing),
      function(.) individual::TargetedEvent$new(parameters$human_population)
    )
    events$pev_epi_doses <- pev_epi_doses
    events$pev_epi_boosters <- pev_epi_boosters
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
  if (!is.null(parameters$mass_pev_timesteps)) {
    events$mass_pev$schedule(parameters$mass_pev_timesteps - 1)
  }
  if (parameters$mda) {
    events$mda_administer$schedule(parameters$mda_timesteps - 1)
  }
  if (parameters$smc) {
    events$smc_administer$schedule(parameters$smc_timesteps - 1)
  }
  if (parameters$tbv) {
    events$tbv_vaccination$schedule(parameters$tbv_timesteps - 1)
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

  # Vaccination event listeners
  if (!is.null(events$mass_pev)) {
    # set up distribution
    events$mass_pev$add_listener(
      create_mass_pev_listener(
        variables,
        events,
        parameters,
        correlations
      )
    )
    attach_pev_dose_listeners(
      variables,
      parameters,
      parameters$mass_pev_timesteps,
      events$mass_pev_doses,
      events$mass_pev_boosters,
      parameters$mass_pev_booster_spacing,
      parameters$mass_pev_booster_coverage,
      parameters$mass_pev_profile_indices,
      'mass',
      renderer
    )
  }

  if (!is.null(events$pev_epi_doses)) {
    attach_pev_dose_listeners(
      variables,
      parameters,
      parameters$pev_epi_timesteps,
      events$pev_epi_doses,
      events$pev_epi_boosters,
      parameters$pev_epi_booster_spacing,
      parameters$pev_epi_booster_coverage,
      parameters$pev_epi_profile_indices,
      'epi',
      renderer
    )
  }

  if (parameters$mda == 1) {
    events$mda_administer$add_listener(create_mda_listeners(
      variables,
      parameters$mda_drug,
      parameters$mda_timesteps,
      parameters$mda_coverages,
      parameters$mda_min_ages,
      parameters$mda_max_ages,
      correlations,
      'mda',
      parameters,
      renderer
    ))
  }

  if (parameters$smc == 1) {
    events$smc_administer$add_listener(create_mda_listeners(
      variables,
      parameters$smc_drug,
      parameters$smc_timesteps,
      parameters$smc_coverages,
      parameters$smc_min_ages,
      parameters$smc_max_ages,
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

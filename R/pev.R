#' @title EPI PEV vaccination process
#'
#' @description schedules individuals to be vaccinated according to the epi
#' strategy
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
create_epi_pev_process <- function(
  variables,
  events,
  parameters,
  correlations,
  coverages,
  timesteps
  ) {
  function(timestep) {
    timestep_index <- match_timestep(ts = timesteps, t = timestep)
    if(timestep_index == 0){
      return()
    }
    coverage <- coverages[timestep_index]
    if(coverage == 0){
      return()
    }
    
    to_vaccinate <- variables$birth$get_index_of(
      set = timestep - parameters$pev_epi_age
    )

    #ignore those who are scheduled for mass vaccination
    if (!is.null(events$mass_pev_doses)) {
      to_vaccinate <- to_vaccinate$and(
        events$mass_pev_doses[[1]]$get_scheduled()$not()
      )
    }

    if (parameters$pev_epi_min_wait == 0) {
      target <- to_vaccinate$to_vector()
    } else {
      not_recently_vaccinated <- variables$last_pev_timestep$get_index_of(
        a = max(timestep - parameters$pev_epi_min_wait, 0),
        b = timestep
      )$not(TRUE)
      target <- to_vaccinate$and(not_recently_vaccinated)$to_vector()
    }

    target <- target[
      sample_intervention(
        target,
       'pev',
        coverage,
        correlations
      )
    ]

    # Update the latest vaccination time
    variables$last_pev_timestep$queue_update(timestep, target)

    schedule_vaccination(
      target,
      events,
      parameters,
      events$pev_epi_doses
    )
  }
}

#' @title mass PEV listener
#'
#' @description schedules individuals to be vaccinated according to the mass
#' strategy
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
create_mass_pev_listener <- function(
  variables,
  events,
  parameters,
  correlations
  ) {
  function(timestep) {
    in_age_group <- individual::Bitset$new(parameters$human_population)
    for (i in seq_along(parameters$mass_pev_min_ages)) {
      min_birth <- timestep - parameters$mass_pev_max_ages[[i]]
      max_birth <- timestep - parameters$mass_pev_min_ages[[i]]
      in_age_group$or(variables$birth$get_index_of(a = min_birth, b = max_birth))
    }
    if (parameters$mass_pev_min_wait == 0) {
      target <- in_age_group
    } else {
      not_recently_vaccinated <- variables$last_pev_timestep$get_index_of(
        a = max(timestep - parameters$mass_pev_min_wait, 0),
        b = timestep
      )$not(TRUE)
      target <- in_age_group$and(not_recently_vaccinated)
    }

    #ignore those who are scheduled for EPI vaccination
    if (!is.null(events$pev_epi_doses)) {
      target <- target$and(
        events$pev_epi_doses[[1]]$get_scheduled()$not()
      )$to_vector()
    } else {
      target <- target$to_vector()
    }
    
    time_index = which(parameters$mass_pev_timesteps == timestep)
    target <- target[
      sample_intervention(
        target,
       'pev',
        parameters$mass_pev_coverages[[time_index]],
        correlations
      )
    ]

    # Update the latest vaccination time
    variables$last_pev_timestep$queue_update(timestep, target)

    # Schedule future doses
    schedule_vaccination(
      target,
      events,
      parameters,
      events$mass_pev_doses
    )
    if (time_index < length(parameters$mass_pev_timesteps)) {
      events$mass_pev$schedule(
        parameters$mass_pev_timesteps[[time_index + 1]] - timestep
      )
    }
  }
}

#' @title Schedule vaccination doses and efficacy
#'
#' @param target vector of individuals to target
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param dose_events a list of dose events to schedule
#' @noRd
schedule_vaccination <- function(
  target,
  events,
  parameters,
  dose_events
  ) {
  if (length(target) > 0) {
    for (d in seq_along(parameters$pev_doses)) {
      dose_events[[d]]$schedule(target, parameters$pev_doses[[d]])
    }
  }
}

#' @title pev efficacy listener
#'
#' @description creates a listener to start pev efficacy in individuals
#'
#' @param variables list of variables in the model
#' @param pev_profile_index the index of the pev profile to introduce
#' @param parameters the model parameters
#' @noRd
create_pev_efficacy_listener <- function(variables, pev_profile_index) {
  function(timestep, target) {
    if (target$size() > 0) {
      variables$last_eff_pev_timestep$queue_update(timestep, target)
      variables$pev_profile$queue_update(pev_profile_index, target)
    }
  }
}

create_pev_booster_listener <- function(
  variables,
  coverage,
  timed_coverage = NULL,
  timed_coverage_timestep = NULL,
  booster_number,
  pev_profile_index,
  next_booster_event,
  next_booster_delay,
  renderer,
  strategy
  ) {
  render_name <- paste0("n_pev_", strategy, "_booster_", booster_number)
  renderer$set_default(render_name, 0)
  force(next_booster_event) # because R lazy evaluation is rubbish
  force(next_booster_delay)
  force(coverage)
  function(timestep, target) {
    if (is.null(timed_coverage)) {
      t_coverage <- 1
    } else {
      t_coverage <- timed_coverage[
        match_timestep(timed_coverage_timestep, timestep)
      ]
    }
    target <- sample_bitset(
      target,
      coverage * t_coverage
    )
    variables$last_pev_timestep$queue_update(timestep, target)
    variables$last_eff_pev_timestep$queue_update(timestep, target)
    variables$pev_profile$queue_update(pev_profile_index, target)
    renderer$render(render_name, target$size(), timestep)

    if (!is.null(next_booster_event)) {
      next_booster_event$schedule(target, next_booster_delay)
    }
  }
}

calculate_pev_antibodies <- function(
  t,
  cs,
  rho,
  ds,
  dl,
  parameters
  ) {
  cs * (
    rho * exp(-t * log(2) / ds) + (
      1 - rho
    ) * exp(-t * log(2) / dl)
  )
}

calculate_pev_efficacy <- function(antibodies, vmax, beta, alpha) {
  vmax * (
    1 - (1 / (
      1 + (antibodies / beta) ** alpha
    ))
  )
}

create_dosage_renderer <- function(renderer, strategy, dose) {
  output_name <- paste0('n_pev_', strategy  ,'_dose_', dose)
  renderer$set_default(output_name, 0)
  function(t, target) renderer$render(output_name, target$size(), t)
}

attach_pev_dose_listeners <- function(
  variables,
  parameters,
  dose_events,
  booster_events,
  booster_delays,
  booster_coverages,
  booster_timed_coverage,
  booster_timed_coverage_timestep,
  pev_profile_indices,
  strategy,
  renderer
  ) {
  # set up dosing
  for (d in seq_along(dose_events)) {
    dose_events[[d]]$add_listener(
      create_dosage_renderer(renderer, strategy, d)
    )
    # update last vaccination on every primary dose
    dose_events[[d]]$add_listener(
      function(t, target) {
        variables$last_pev_timestep$queue_update(t, target)
      }
    )
    if (d == length(dose_events)) {
      dose_events[[d]]$add_listener(
        create_pev_efficacy_listener(
          variables,
          pev_profile_indices[[1]]
        )
      )
      if (length(booster_events) > 0) {
        seasonal_boosters <- FALSE
        if (!is.null(parameters$pev_epi_seasonal_boosters)) {
          seasonal_boosters <- parameters$pev_epi_seasonal_boosters
        }
        if (seasonal_boosters) {
          dose_events[[d]]$add_listener(
            create_seasonal_booster_scheduler(
              booster_events[[1]],
              booster_delays[[1]],
              parameters
            )
          )
        } else  {
          dose_events[[d]]$add_listener(
            individual::reschedule_listener(
              booster_events[[1]],
              booster_delays[[1]]
            )
          )
        }
      }
    }
  }

  # set up boosters
  for (b in seq_along(booster_events)) {
    if (b == length(booster_events)) {
      next_booster_event <- NULL
      next_booster_delay <- NULL
    } else {
      next_booster_event <- booster_events[[b + 1]]
      next_booster_delay <- diff(
        booster_delays[c(b, b + 1)]
      )
    }

    booster_events[[b]]$add_listener(
      create_pev_booster_listener(
        variables = variables,
        coverage = booster_coverages[[b]],
        timed_coverage = booster_timed_coverage,
        timed_coverage_timestep = booster_timed_coverage_timestep,
        booster_number = b,
        pev_profile_index = pev_profile_indices[[b + 1]],
        next_booster_event = next_booster_event,
        next_booster_delay = next_booster_delay,
        renderer = renderer,
        strategy = strategy
      )
    )
  }
}

create_seasonal_booster_scheduler <- function(
  booster_event,
  booster_delay,
  parameters
  ) {
  function(timestep, target) {
    delay <- booster_delay - timestep %% 365
    if (delay < 0) {
      delay <- delay + 365
    }
    if (delay <= parameters$pev_epi_min_wait) {
      delay <- delay + 365
    }
    booster_event$schedule(target, delay)
  }
}

sample_pev_param <- function(profile_index, profile_list, param_name) {
  mu <- vnapply(profile_list, function(p) p[[param_name]][[1]])
  sigma <- vnapply(profile_list, function(p) p[[param_name]][[2]])
  rnorm(length(profile_index), mu[profile_index], sigma[profile_index])
}

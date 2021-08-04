#' @title RTS,S EPI vaccination process
#'
#' @description schedules individuals to be vaccinated according to the epi
#' strategy
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
create_rtss_epi_process <- function(
  variables,
  events,
  parameters,
  correlations
  ) {
  function(timestep) {
    if (timestep <- parameters$rtss_epi_timestep) {
      return()
    }
    to_vaccinate <- variables$birth$get_index_of(
      set = timestep - parameters$rtss_epi_age
    )
    not_recently_vaccinated <- variables$rtss_vaccinated(
      a = timestep - parameters$rtss_min_wait,
      b = timestep
    )$not()
    target <- to_vaccinate$and(not_recently_vaccinated)$to_vector()
    target <- target[
      sample_intervention(
        target,
       'rtss',
        parameters$rtss_epi_coverage,
        correlations
      )
    ]
    if (target$size() > 0) {
      schedule_vaccination(
        target,
        variables,
        events,
        parameters
      )
    }
  }
}

#' @title RTS,S mass vaccination listener
#'
#' @description schedules individuals to be vaccinated according to the mass
#' strategy
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
create_rtss_mass_listener <- function(
  variables,
  events,
  parameters,
  correlations
  ) {
  function(timestep) {
    time_index = which(parameters$rtss_timesteps == timestep)
    not_recently_vaccinated <- variables$rtss_vaccinated(
      a = timestep - parameters$rtss_min_wait,
      b = timestep
    )$not()
    in_age_group <- individual::Bitset$new(parameters$human_population)
    for (i in seq_along(parameters$rtss_min_ages)) {
      min_birth <- timestep - parameters$rtss_max_ages[[i]]
      max_birth <- timestep - parameters$rtss_min_ages[[i]]
      in_age_group$or(variables$birth$get_index_of(a = min_birth, b = max_birth))
    }
    target <- in_age_group$and(not_recently_vaccinated)$to_vector()
    target <- target[
      sample_intervention(
        target,
       'rtss',
        parameters$rtss_coverages[[time_index]],
        correlations
      )
    ]
    schedule_vaccination(target, variables, events, parameters)
    if (time_index < length(parameters$rtss_timesteps)) {
      events$rtss_mass_vaccination$schedule(
        parameters$rtss_timesteps[[time_index + 1]] - timestep
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
#' @noRd
schedule_vaccination <- function(
  target,
  events,
  parameters
  ) {
  if (target$size() > 0) {
    for (d in seq_along(parameters$rtss_doses)) {
      events$rtss_doses[[d]]$schedule(target, parameters$rtss_doses[[d]])
    }
  }
}

#' @title RTS,S efficacy listener
#'
#' @description creates a listener to start vaccine efficacy in individuals
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @noRd
create_rtss_efficacy_listener <- function(
  variables,
  events,
  parameters,
  boosters,
  booster_coverage
  ) {
  function(timestep, target) {
    if (target$size() > 0) {
      variables$rtss_vaccinated$queue_update(
        timestep,
        target
      )
      if (length(parameters$rtss_boosters) > 0) {
        events$rtss_booster$schedule(
          target$sample(booster_coverage[[1]]),
          boosters[[1]]
        )
      }
    }
  }
}

create_rtss_booster_listener <- function(
  variables,
  events,
  parameters,
  boosters,
  booster_coverage
  ) {
  function(timestep, target) {
    variables$rtss_cs$queue_update(
      exp(
        parameters$rtss_cs_boost[[1]] + parameters$rtss_cs_boost[[2]] * rnorm(target$size())
      ),
      target
    )

    variables$rtss_rho$queue_update(
      invlogit(
        parameters$rtss_rho_boost[[1]] + parameters$rtss_rho_boost[[2]] * rnorm(target$size())
      ),
      target
    )

    variables$rtss_boosted$queue_update(timestep, target)

    vaccinated <- variables$rtss_vaccinated$get_values(target)
    for (v in unique(vaccinated)) {
      for (i in seq_len(length(boosters) - 1)) {
        to_boost <- which(
          vaccinated + boosters[[i]] == timestep
        )

        if (length(to_boost) > 0) {
          to_boost <- to_boost[bernoulli(
            length(to_boost),
            booster_coverage[[i + 1]]
          )]
          events$rtss_booster$schedule(
            to_boost,
            v + boosters[[i + 1]] - timestep
          )
        }
      }
    }
  }
}

calculate_rtss_antibodies <- function(
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

calculate_rtss_efficacy <- function(antibodies, parameters) {
  parameters$rtss_vmax * (
    1 - (1 / (
      1 + (antibodies / parameters$rtss_beta) ** parameters$rtss_alpha
    ))
  )
}

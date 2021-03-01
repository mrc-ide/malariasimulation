#' @title RTS,S vaccination
#' @description creates a listener which models vaccination according to the
#' strategy from `set_rtss` and correlation parameters from
#' `get_correlation_parameters`
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @param renderer rendering object
#' @noRd
create_rtss_vaccination_listener <- function(
  variables,
  events,
  parameters,
  correlations,
  renderer
  ) {
  function(timestep) {
    time_index = which(parameters$rtss_timesteps == timestep)
    age <- get_age(variables$birth$get_values(), timestep)
    not_vaccinated <- variables$rtss_vaccinated$get_values() == -1
    target_indices <- rep(FALSE, length(age))
    for (i in seq_along(parameters$rtss_min_ages)) {
      target_indices <- target_indices | (
        age >= parameters$rtss_min_ages[[i]] &
          age <= parameters$rtss_max_ages[[i]]
      )
    }
    target <- which(target_indices & not_vaccinated)
    target <- target[
      sample_intervention(
        target,
        'rtss',
        parameters$rtss_coverages[[time_index]],
        correlations
      )
    ]
    renderer$render('n_vaccinated', length(target), timestep)
    if (length(target) > 0) {
      variables$rtss_vaccinated$queue_update(
        timestep,
        target
      )
    }
    if (time_index < length(parameters$rtss_timesteps)) {
      events$rtss_vaccination$schedule(
        parameters$rtss_timesteps[[time_index + 1]] - timestep
      )
    }
    if (length(parameters$rtss_boosters) > 0) {
      events$rtss_booster$schedule(
        target[bernoulli(length(target), parameters$rtss_booster_coverage[[1]])],
        parameters$rtss_boosters[[1]]
      )
    }
  }
}

create_rtss_booster_listener <- function(variables, events, parameters) {
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

    variables$rtss_boosted$queue_update(
      timestep,
			target
		)

    vaccinated <- variables$rtss_vaccinated$get_values(target)
    for (v in unique(vaccinated)) {
      for (i in seq_len(length(parameters$rtss_boosters) - 1)) {
        to_boost <- which(
          vaccinated + parameters$rtss_boosters[[i]] == timestep
        )

        if (length(to_boost) > 0) {
          to_boost <- to_boost[bernoulli(
            length(to_boost),
            parameters$rtss_booster_coverage[[i + 1]]
          )]
          events$rtss_booster$schedule(
            to_boost,
            v + parameters$rtss_boosters[[i + 1]] - timestep
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

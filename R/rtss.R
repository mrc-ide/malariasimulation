#' @title RTS,S vaccination
#' @description creates a listener which models vaccination according to the
#' strategy from `set_rtss` and correlation parameters from
#' `get_correlation_parameters`
#'
#' @param human the handle for the human individual
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
create_rtss_vaccination_listener <- function(
  human,
  variables,
  events,
  parameters,
  correlations
  ) {
  function(api, target) {
    timestep <- api$get_timestep()
    age <- get_age(
      api$get_variable(human, variables$birth),
      timestep
    )
    not_vaccinated <- api$get_variable(human, variables$rtss_vaccinated) == -1
    target_indices <- rep(FALSE, length(age))
    for (i in seq_along(parameters$rtss_min_ages)) {
      target_indices <- target_indices | (
        age >= parameters$rtss_min_ages[[i]] &
          age <= parameters$rtss_max_ages[[i]]
      )
    }
    target <- which(target_indices & not_vaccinated)
    target <- target[
      sample_intervention(target, 'rtss', parameters$rtss_coverage, correlations)
    ]
    api$render('n_vaccinated', length(target))
    if (length(target) > 0) {
      api$queue_variable_update(
        human,
        variables$rtss_vaccinated,
        timestep,
        target
      )
    }
    if (timestep + parameters$rtss_frequency <= parameters$rtss_end) {
      api$schedule(events$rtss_vaccination, c(1), parameters$rtss_frequency)
    }
    if (length(parameters$rtss_boosters) > 0) {
      api$schedule(
        events$rtss_booster,
        target[bernoulli(length(target), parameters$rtss_booster_coverage[[1]])],
        parameters$rtss_boosters[[1]]
      )
    }
  }
}

create_rtss_booster_listener <- function(human, variables, events, parameters) {
  function(api, target) {
		api$queue_variable_update(
			human,
			variables$rtss_cs,
			exp(
        parameters$rtss_cs_boost[[1]] + parameters$rtss_cs_boost[[2]] * rnorm(length(target))
      ),
      target
		)

		api$queue_variable_update(
			human,
			variables$rtss_rho,
			invlogit(
        parameters$rtss_rho_boost[[1]] + parameters$rtss_rho_boost[[2]] * rnorm(length(target))
      ),
			target
		)

    timestep <- api$get_timestep()
    api$queue_variable_update(
			human,
			variables$rtss_boosted,
      timestep,
			target
		)

    vaccinated <- api$get_variable(human, variables$rtss_vaccinated, target)
    for (v in unique(vaccinated)) {
      for (i in seq_len(length(parameters$rtss_boosters) - 1)) {
        to_boost <- which(
          vaccinated + parameters$rtss_boosters[[i]] == timestep
        )

        if (length(to_boost) > 0) {
          to_boost <- to_boost[
            bernoulli(length(to_boost), parameters$rtss_booster_coverage[[i + 1]])
          ]
          api$schedule(
            events$rtss_booster,
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

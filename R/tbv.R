#' @title Calculate the effect of TBV on infectivity
#' @description Returns a vector of human infectivity towards mosquitoes
#' accounting for the reduction in transmission due to vaccination
#'
#' @param timestep current timestep
#' @param infectivity a vector of raw infectivities
#' @param variables the available variables
#' @param parameters model parameters
#' @noRd
account_for_tbv <- function(
  timestep,
  infectivity,
  variables,
  parameters
) {
  time_vaccinated <- variables$tbv_vaccinated$get_values()
  vaccinated <- which(time_vaccinated != -1)
  affected_states <- c('U', 'A', 'D', 'Tr')
  mx <- parameters[c('tbv_mu', 'tbv_ma', 'tbv_md', 'tbv_mt')]
  for (i in seq_along(affected_states)) {
    in_state <- variables$state$get_index_of(affected_states[[i]])$to_vector()
    vaccinated_in_state <- intersect(vaccinated, in_state)
    antibodies <- calculate_tbv_antibodies(
      timestep - time_vaccinated[vaccinated_in_state]#,
      #parameters$tbv_tau,
      #parameters$tbv_rho,
      #parameters$tbv_ds,
      #parameters$tbv_dl
    )
    tra <- calculate_TRA(
      #parameters$tbv_tra_mu,
      #parameters$tbv_gamma1,
      #parameters$tbv_gamma2,
      parameters$tbv_hill,
      parameters$tbv_EC50,
      antibodies
    )
    tba <- calculate_TBA(
      mx[[i]],
      parameters$tbv_k,
      tra
    )
    infectivity[vaccinated_in_state] <- infectivity[vaccinated_in_state] * (
      1 - tba
    )
  }
  infectivity
}


#' @title Distribute TBV vaccine
#'
#' @param timestep current timestep
#' @param variables the available variables
#' @param events the available events
#' @param parameters model parameters
#' @param correlations model correlations
#' @param renderer object for model outputs
#' @noRd
create_tbv_listener <- function(variables, events, parameters, correlations, renderer) {
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
}

calculate_tbv_antibodies <- function(t) { #this'll need updating. Add more arguments??
  #tau * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl))
  199*(0.36*exp(-t/6) + 0.44*exp(-t/37) + 0.2*exp(-t/50))
}

calculate_TRA <- function(hill, EC50, antibodies) { #this'll need updating
  #numerator <- (antibodies / mu)^gamma1
  #numerator / (numerator + gamma2)
  (antibodies^hill)/((antibodies^hill) + (EC50^hill))
}

calculate_TBA <- function(mx, k, tra) { #same
  offset <- (k / (k + mx)) ^ k;
  scale <- 1. / (1. - offset);
  tra_transformation <- (k / (k + mx * (1 - tra))) ^ k;
  scale * (tra_transformation - offset)
}

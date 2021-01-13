#' @title Calculate overall infections for bitten humans
#' @description
#' Sample infected humans given prophylaxis and vaccination
#' @param variables a list of all of the model variables
#' @param bitten_humans indices of bitten humans
#' @param ib vector of pre-erythrocytic immunity levels
#' @param timestep current timestep
#' @importFrom stats dweibull
calculate_infections <- function(
  variables,
  bitten_humans,
  ib,
  parameters,
  timestep
  ) {
  source_humans <- variables$state$get_index_of(c('S', 'A', 'U'))
  bitten_bitset <- individual::Bitset$new(source_humans$max_size)
  bitten_bitset$insert(bitten_humans)
  source_humans <- source_humans$and(bitten_bitset)$to_vector()

  b <- blood_immunity(ib[source_humans], parameters)

  # calculate prophylaxis
  prophylaxis <- rep(0, length(source_humans))
  drug <- variables$drug$get_values(source_humans)
  medicated <- (drug > 0)
  if (any(medicated)) {
    drug <- drug[medicated]
    drug_time <- variables$drug_time$get_values(source_humans[medicated])
    prophylaxis[medicated] <- dweibull(
      timestep - drug_time,
      parameters$drug_prophylaxis_shape[drug],
      parameters$drug_prophylaxis_scale[drug]
    )
  }

  # calculate vaccine efficacy
  vaccine_efficacy <- rep(0, length(source_humans))
  vaccine_times <- pmax(
    variables$rtss_vaccinated$get_values(source_humans),
    variables$rtss_boosted$get_values(source_humans)
  )
  vaccinated <- which(vaccine_times > -1)
  vaccinated_index <- source_humans[vaccine_times > -1]
  antibodies <- calculate_rtss_antibodies(
    timestep - vaccine_times[vaccinated],
    variables$rtss_cs$get_values(vaccinated_index),
    variables$rtss_rho$get_values(vaccinated_index),
    variables$rtss_ds$get_values(vaccinated_index),
    variables$rtss_dl$get_values(vaccinated_index),
    parameters
  )
  vaccine_efficacy[vaccinated] <- calculate_rtss_efficacy(antibodies, parameters)

  source_humans[bernoulli_multi_p(
    b * (1 - prophylaxis) * (1 - vaccine_efficacy)
  )]
}

#' @title Calculate clinical infections
#' @description
#' Sample clinical infections from all infections
#' @param variables a list of all of the model variables
#' @param infections index of infected humans
#' @param parameters model parameters
calculate_clinical_infections <- function(variables, infections, parameters, timestep) {
  ica <- variables$ica$get_values(infections)
  icm <- variables$icm$get_values(infections)

  if (length(infections) > 0) {
    boost_immunity(
      variables$ica,
      infections,
      ica,
      variables$last_boosted_ica,
      timestep,
      parameters$uc
    )
    boost_immunity(
      variables$id,
      infections,
      variables$id$get_values(infections),
      variables$last_boosted_id,
      timestep,
      parameters$ud
    )
  }

  phi <- clinical_immunity(ica, icm, parameters)
  infections[bernoulli_multi_p(phi)]
}

#' @title Calculate severe infections
#' @description
#' Sample severely infected humans from clinically infected
#' @param api simulation api
#' @param clinical_infections indices of clinically infected humans
#' @param infection_age ages of individuals in `clinical_infections` (timesteps)
#' @param human handle for humans
#' @param variables a list of all of the model variables
#' @param infections indices of all infected humans (for immunity boosting)
update_severe_disease <- function(
  api,
  clinical_infections,
  infection_age,
  human,
  variables,
  infections
  ) {
  if (length(clinical_infections) > 0) {
    parameters <- api$get_parameters()
    iva <- api$get_variable(human, variables$iva, clinical_infections)
    theta <- severe_immunity(
      infection_age,
      iva,
      api$get_variable(human, variables$ivm, clinical_infections),
      parameters
    )
    develop_severe <- bernoulli_multi_p(theta)
    api$queue_variable_update(
      human,
      variables$is_severe,
      develop_severe,
      clinical_infections
    )
    boost_immunity(
      api,
      human,
      variables$iva,
      infections,
      iva,
      variables$last_boosted_iva,
      api$get_timestep(),
      parameters$uv
    )
  }
}

#' @title Calculate treated humans
#' @description
#' Sample treated humans from the clinically infected
#' @param variables a list of all of the model variables
#' @param clinical_infections indices of clinically infected humans
#' @param recovery the recovery event
#' @param parameters model parameters
#' @param timestep
calculate_treated <- function(
  variables,
  clinical_infections,
  recovery,
  parameters,
  timestep
  ) {
  if (length(parameters$clinical_treatment_coverages) == 0) {
    return(numeric(0))
  }

  seek_treatment <- bernoulli(length(clinical_infections), parameters$ft)
  n_treat <- length(seek_treatment)
  drugs <- parameters$clinical_treatment_drugs[
    sample.int(
      length(parameters$clinical_treatment_drugs),
      n_treat,
      prob = parameters$clinical_treatment_coverages,
      replace = TRUE
    )
  ]

  successful <- bernoulli_multi_p(parameters$drug_efficacy[drugs])
  treated_index <- clinical_infections[seek_treatment][successful]

  # Update those who have been treated
  if (length(treated_index) > 0) {
    variables$state$queue_update('Tr', treated_index)
    variables$infectivity$queue_update(
      parameters$cd * parameters$drug_rel_c[drugs[successful]],
      treated_index
    )
    variables$drug$queue_update(
      drugs[successful],
      treated_index
    )
    variables$drug_time$queue_update(
      timestep,
      treated_index
    )
    recovery$schedule(
      treated_index,
      log_uniform(length(treated_index), parameters$dt)
    )
  }
  treated_index
}

#' @title Schedule infections
#' @description
#' Schedule infections in humans after the incubation period
#' @param events a list of all of the model events
#' @param clinical_infections indices of clinically infected humans
#' @param treated indices of treated humans
#' @param infections indices of infected humans
schedule_infections <- function(
  events,
  clinical_infections,
  treated,
  infections,
  parameters
  ) {
  scheduled_for_infection <- events$infection$get_scheduled()
  treated_b <- individual::Bitset$new(scheduled_for_infection$max_size)
  treated_b$insert(treated)
  excluded <- scheduled_for_infection$and(treated)

  infections
  all_new_infections <- setdiff(infections, excluded)
  to_infect <- all_new_infections %in% clinical_infections
  to_infect_asym <- !to_infect

  infection_times <- log_uniform(length(all_new_infections), parameters$de)

  if(sum(to_infect) > 0) {
    api$schedule(
      events$clinical_infection,
      all_new_infections[to_infect],
      infection_times[to_infect]
    )
  }

  if(sum(to_infect_asym) > 0) {
    api$schedule(
      events$asymptomatic_infection,
      all_new_infections[to_infect_asym],
      infection_times[to_infect_asym]
    )
  }

  if(length(all_new_infections) > 0) {
    api$schedule(events$infection, all_new_infections, infection_times)
    api$clear_schedule(events$subpatent_infection, all_new_infections)
    api$clear_schedule(events$recovery, all_new_infections)
  }
}

# =================
# Utility functions
# =================
boost_immunity <- function(
  api,
  human,
  immunity_variable,
  exposed_index,
  exposed_values,
  last_boosted_variable,
  timestep,
  delay
  ) {
  # record who can be boosted
  last_boosted <- api$get_variable(human, last_boosted_variable, exposed_index)
  to_boost <- (timestep - last_boosted) >= delay | (last_boosted == -1)
  if (sum(to_boost) > 0) {
    # boost the variable
    api$queue_variable_update(
      human,
      immunity_variable,
      exposed_values[to_boost] + 1,
      exposed_index[to_boost]
    )
    # record last boosted
    api$queue_variable_update(
      human,
      last_boosted_variable,
      timestep,
      exposed_index[to_boost]
    )
  }
}

# Implemented from Winskill 2017 - Supplementary Information page 4
clinical_immunity <- function(acquired_immunity, maternal_immunity, parameters) {
  parameters$phi0 * (
    parameters$phi1 +
      (1 - parameters$phi1) /
      (
        1 + (
          (acquired_immunity + maternal_immunity) / parameters$ic0
        ) ** parameters$kc
      )
  )
}

# Implemented from Winskill 2017 - Supplementary Information page 5
severe_immunity <- function(age, acquired_immunity, maternal_immunity, parameters) {
  fv <- 1 - (1 - parameters$fv0) / (
    1 + (age / parameters$av) ** parameters$gammav
  )
  parameters$theta0 * (parameters$theta1 + (1 - parameters$theta1) / (
    1 + fv * (
      (acquired_immunity + maternal_immunity) / parameters$iv0) ** parameters$kv
    )
  )
}

# Implemented from Winskill 2017 - Supplementary Information page 5
# NOTE: I believe there is a typo on equation (9) the + 1 in
# the denominator is in the wrong place
# NOTE: Also dmin = d1, according to the equilibrium solution `main.R:132`
probability_of_detection <- function(age, immunity, parameters) {
  fd <- 1 - (1 - parameters$fd0) / (
    1 + (age / parameters$ad) ** parameters$gammad
  )
  parameters$d1 + (1 - parameters$d1) / (
    1 + fd * (immunity / parameters$id0) ** parameters$kd
  )
}

# Implemented from Winskill 2017 - Supplementary Information page 6
# NOTE: there appears to be a typo on line 114, should be (cd - cu)
asymptomatic_infectivity <- function(age, immunity, parameters) {
  q <- probability_of_detection(age, immunity, parameters)
  parameters$cu + (parameters$cd - parameters$cu) * q ** parameters$gamma1
}

# Unique biting rate (psi) for a human of a given age
unique_biting_rate <- function(age, parameters) {
  1 - parameters$rho * exp(- age / parameters$a0)
}

# Implemented from Winskill 2017 - Supplementary Information page 4
blood_immunity <- function(ib, parameters) {
  parameters$b0 * (
    parameters$b1 +
      (1 - parameters$b1) /
      (1 + (ib / parameters$ib0) ** parameters$kb)
  )
}

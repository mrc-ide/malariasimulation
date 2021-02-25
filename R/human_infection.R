#' @title Calculate overall infections for bitten humans
#' @description
#' Sample infected humans given prophylaxis and vaccination
#' @param variables a list of all of the model variables
#' @param bitten_humans bitset of bitten humans
#' @param ib vector of pre-erythrocytic immunity levels
#' @param timestep current timestep
#' @importFrom stats dweibull
calculate_infections <- function(
  variables,
  bitten_humans,
  parameters,
  timestep
  ) {
  source_humans <- variables$state$get_index_of(
    c('S', 'A', 'U'))$and(bitten_humans)

  b <- blood_immunity(variables$ib$get_values(source_humans), parameters)

  source_vector <- source_humans$to_vector()

  # calculate prophylaxis
  prophylaxis <- rep(0, length(source_vector))
  drug <- variables$drug$get_values(source_vector)
  medicated <- (drug > 0)
  if (any(medicated)) {
    drug <- drug[medicated]
    drug_time <- variables$drug_time$get_values(source_vector[medicated])
    prophylaxis[medicated] <- dweibull(
      timestep - drug_time,
      parameters$drug_prophylaxis_shape[drug],
      parameters$drug_prophylaxis_scale[drug]
    )
  }

  # calculate vaccine efficacy
  vaccine_efficacy <- rep(0, length(source_vector))
  vaccine_times <- pmax(
    variables$rtss_vaccinated$get_values(source_vector),
    variables$rtss_boosted$get_values(source_vector)
  )
  vaccinated <- which(vaccine_times > -1)
  if (length(vaccinated) > 0) {
    vaccinated_index <- source_vector[vaccine_times > -1]
    antibodies <- calculate_rtss_antibodies(
      timestep - vaccine_times[vaccinated],
      variables$rtss_cs$get_values(vaccinated_index),
      variables$rtss_rho$get_values(vaccinated_index),
      variables$rtss_ds$get_values(vaccinated_index),
      variables$rtss_dl$get_values(vaccinated_index),
      parameters
    )
    vaccine_efficacy[vaccinated] <- calculate_rtss_efficacy(antibodies, parameters)
  }

  bitset_at(
    source_humans,
    bernoulli_multi_p(b * (1 - prophylaxis) * (1 - vaccine_efficacy))
  )
}

#' @title Calculate clinical infections
#' @description
#' Sample clinical infections from all infections
#' @param variables a list of all of the model variables
#' @param infections bitset of infected humans
#' @param parameters model parameters
#' @param timestep current timestep
calculate_clinical_infections <- function(variables, infections, parameters, timestep) {
  ica <- variables$ica$get_values(infections)
  icm <- variables$icm$get_values(infections)

  if (infections$size() > 0) {
    boost_immunity(
      variables$ica,
      infections,
      variables$last_boosted_ica,
      timestep,
      parameters$uc
    )
    boost_immunity(
      variables$id,
      infections,
      variables$last_boosted_id,
      timestep,
      parameters$ud
    )
  }

  phi <- clinical_immunity(ica, icm, parameters)
  bitset_at(infections, bernoulli_multi_p(phi))
}

#' @title Calculate severe infections
#' @description
#' Sample severely infected humans from clinically infected
#' @param timestep current timestep
#' @param clinical_infections indices of clinically infected humans
#' @param variables a list of all of the model variables
#' @param infections indices of all infected humans (for immunity boosting)
#' @param parameters model parameters
update_severe_disease <- function(
  timestep,
  clinical_infections,
  variables,
  infections,
  parameters
  ) {
  if (clinical_infections$size() > 0) {
    age <- get_age(variables$birth$get_values(clinical_infections), timestep)
    iva <- variables$iva$get_values(clinical_infections)
    theta <- severe_immunity(
      age,
      iva,
      variables$ivm$get_values(clinical_infections),
      parameters
    )
    develop_severe <- bernoulli_multi_p(theta)
    variables$is_severe$queue_update(
      'yes',
      bitset_at(clinical_infections, develop_severe)
    )
    boost_immunity(
      variables$iva,
      infections,
      variables$last_boosted_iva,
      timestep,
      parameters$uv
    )
  }
}

#' @title Calculate treated humans
#' @description
#' Sample treated humans from the clinically infected
#' @param variables a list of all of the model variables
#' @param clinical_infections a bitset of clinically infected humans
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
    return(individual::Bitset$new(parameters$human_population))
  }

  seek_treatment <- sample_bitset(clinical_infections, parameters$ft)
  n_treat <- seek_treatment$size()
  drugs <- parameters$clinical_treatment_drugs[
    sample.int(
      length(parameters$clinical_treatment_drugs),
      n_treat,
      prob = parameters$clinical_treatment_coverages,
      replace = TRUE
    )
  ]

  successful <- bernoulli_multi_p(parameters$drug_efficacy[drugs])
  treated_index <- bitset_at(seek_treatment, successful)

  # Update those who have been treated
  if (treated_index$size() > 0) {
    successful_vector <- successful$to_vector()
    variables$state$queue_update('Tr', treated_index)
    variables$infectivity$queue_update(
      parameters$cd * parameters$drug_rel_c[drugs[successful_vector]],
      treated_index
    )
    variables$drug$queue_update(
      drugs[successful_vector],
      treated_index
    )
    variables$drug_time$queue_update(
      timestep,
      treated_index
    )
    recovery$schedule(
      treated_index,
      log_uniform(treated_index$size(), parameters$dt)
    )
  }
  treated_index
}

#' @title Schedule infections
#' @description
#' Schedule infections in humans after the incubation period
#' @param events a list of all of the model events
#' @param clinical_infections bitset of clinically infected humans
#' @param treated bitset of treated humans
#' @param infections bitset of infected humans
schedule_infections <- function(
  events,
  clinical_infections,
  treated,
  infections,
  parameters
  ) {
  included <- events$infection$get_scheduled()$or(treated)$not()

  to_infect <- clinical_infections$and(included)
  to_infect_asym <- clinical_infections$not()$and(infections)$and(included)

  # change to symptomatic
  if(to_infect$size() > 0) {
    infection_times <- log_uniform(to_infect$size(), parameters$de)
    events$clinical_infection$schedule(
      to_infect,
      infection_times
    )
    events$infection$schedule(to_infect, infection_times)
    events$subpatent_infection$clear_schedule(to_infect)
    events$recovery$clear_schedule(to_infect)
  }

  if(to_infect_asym$size() > 0) {
    infection_times <- log_uniform(to_infect_asym$size(), parameters$de)
    events$asymptomatic_infection$schedule(
      to_infect_asym,
      infection_times
    )
    events$infection$schedule(to_infect_asym, infection_times)
    events$subpatent_infection$clear_schedule(to_infect_asym)
    events$recovery$clear_schedule(to_infect_asym)
  }
}

# =================
# Utility functions
# =================
boost_immunity <- function(
  immunity_variable,
  exposed_index,
  last_boosted_variable,
  timestep,
  delay
  ) {
  # record who can be boosted
  exposed_index_vector <- exposed_index$to_vector()
  last_boosted <- last_boosted_variable$get_values(exposed_index)
  to_boost <- (timestep - last_boosted) >= delay | (last_boosted == -1)
  exposed_to_boost <- exposed_index_vector[to_boost]
  if (sum(to_boost) > 0) {
    # boost the variable
    immunity_variable$queue_update(
      immunity_variable$get_values(exposed_to_boost) + 1,
      exposed_to_boost
    )
    # record last boosted
    last_boosted_variable$queue_update(
      timestep,
      exposed_to_boost
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

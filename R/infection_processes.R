#' @title Infection process
#' @description
#' This is the process of infection for humans. It results in human future state
#' changes for infected humans and boosts in immunity.
#' @param individuals a list of individuals in the model
#' @param states a list of all of the model states
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
create_infection_process <- function(
  individuals,
  states,
  variables,
  events
  ) {
  function(api) {
    human <- individuals$human
    parameters <- api$get_parameters()
    timestep <- api$get_timestep()

    # Calculate EIR
    age <- get_age(api$get_variable(human, variables$birth), api$get_timestep())
    epsilon <- eir_from_api(api, individuals, states, variables, age)

    api$render("mean_EIR", mean(epsilon))

    bitten_humans <- which(bernoulli_multi_p(length(epsilon), epsilon))
    ib <- api$get_variable(human, variables$ib)
    if (length(bitten_humans) > 0) {
      boost_immunity(
        api,
        human,
        variables$ib,
        bitten_humans,
        ib[bitten_humans],
        variables$last_boosted_ib,
        timestep,
        parameters$ub
      )
    }

    # Calculate Infected
    infected_humans <- calculate_infections(
      api,
      human,
      states,
      variables,
      bitten_humans,
      ib
    )

    clinical_infections <- calculate_clinical_infections(
      api,
      human,
      variables,
      infected_humans
    )

    if (parameters$severe_enabled) {
      update_severe_disease(
        api,
        clinical_infections,
        age[clinical_infections],
        human,
        variables,
        infected_humans
      )
    }

    treated <- calculate_treated(
      api,
      human,
      states,
      variables,
      clinical_infections
    )

    schedule_infections(
      api,
      events,
      clinical_infections,
      treated,
      infected_humans
    )
  }
}

eir_from_api <- function(api, individuals, states, variables, age) {
  parameters <- api$get_parameters()
  source_mosquitos <- api$get_state(individuals$mosquito, states$Im)
  infectivity <- vector_infectivity(
    api$get_variable(
      individuals$mosquito,
      variables$mosquito_variety,
      source_mosquitos
    ),
    parameters
  )

  eir(
    age,
    api$get_variable(individuals$human, variables$zeta),
    infectivity,
    parameters
  )
}

#' @importFrom stats dweibull
calculate_infections <- function(
  api,
  human,
  states,
  variables,
  bitten_humans,
  ib
  ) {
  parameters <- api$get_parameters()
  source_humans <- intersect(
    api$get_state(human, states$S, states$A, states$U),
    bitten_humans
  )
  b <- blood_immunity(ib[source_humans], parameters)

  # calculate prophylaxis
  prophylaxis <- rep(0, length(source_humans))
  drug <- api$get_variable(human, variables$drug, source_humans)
  medicated <- (drug > 0)
  if (any(medicated)) {
    drug <- drug[medicated]
    drug_time <- api$get_variable(
      human,
      variables$drug_time,
      source_humans[medicated]
    )
    prophylaxis[medicated] <- dweibull(
      api$get_timestep() - drug_time,
      parameters$drug_prophylaxis_shape[drug],
      parameters$drug_prophylaxis_scale[drug]
    )
  }

  # calculate vaccine efficacy
  vaccine_efficacy <- rep(0, length(source_humans))
  vaccine_times <- pmax(
    api$get_variable(human, variables$rtss_vaccinated, source_humans),
    api$get_variable(human, variables$rtss_boosted, source_humans)
  )
	vaccinated <- which(vaccine_times > -1)
  vaccinated_index <- source_humans[vaccine_times > -1]
  antibodies <- calculate_rtss_antibodies(
    api$get_timestep() - vaccine_times[vaccinated],
    api$get_variable(human, variables$rtss_cs, vaccinated_index),
    api$get_variable(human, variables$rtss_rho, vaccinated_index),
    api$get_variable(human, variables$rtss_ds, vaccinated_index),
    api$get_variable(human, variables$rtss_dl, vaccinated_index),
    parameters
  )
  vaccine_efficacy[vaccinated] <- calculate_rtss_efficacy(antibodies, parameters)

  source_humans[bernoulli_multi_p(
    length(source_humans),
    b * (1 - prophylaxis) * (1 - vaccine_efficacy)
  )]
}

calculate_clinical_infections <- function(api, human, variables, infections) {
  ica <- api$get_variable(human, variables$ica, infections)
  icm <- api$get_variable(human, variables$icm, infections)
  parameters <- api$get_parameters()

  if (length(infections) > 0) {
    timestep <- api$get_timestep()
    boost_immunity(
      api,
      human,
      variables$ica,
      infections,
      ica,
      variables$last_boosted_ica,
      timestep,
      parameters$uc
    )
    boost_immunity(
      api,
      human,
      variables$id,
      infections,
      api$get_variable(human, variables$id, infections),
      variables$last_boosted_id,
      timestep,
      parameters$ud
    )
  }

  phi <- clinical_immunity(ica, icm, parameters)
  infections[bernoulli_multi_p(length(infections), phi)]
}

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
    iva <- api$get_variable(human, variables$iva, infections)
    theta <- severe_immunity(
      infection_age,
      iva,
      api$get_variable(human, variables$ivm, clinical_infections),
      parameters
    )
    develop_severe <- bernoulli_multi_p(length(clinical_infections), theta)
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

calculate_treated <- function(
  api,
  human,
  states,
  variables,
  clinical_infections
  ) {
  parameters <- api$get_parameters()
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

  successful <- bernoulli_multi_p(n_treat, parameters$drug_efficacy[drugs])
  treated_index <- clinical_infections[seek_treatment][successful]

  # Update those who have been treated
  if (length(treated_index) > 0) {
    api$queue_state_update(human, states$Tr, treated_index)
    api$queue_variable_update(
      human,
      variables$infectivity,
      parameters$cd * parameters$drug_rel_c[drugs[successful]],
      treated_index
    )
    api$queue_variable_update(
      human,
      variables$drug,
      drugs[successful],
      treated_index
    )
    api$queue_variable_update(
      human,
      variables$drug_time,
      api$get_timestep(),
      treated_index
    )
  }
  treated_index
}

schedule_infections <- function(
  api,
  events,
  clinical_infections,
  treated,
  infections
  ) {
  parameters <- api$get_parameters()
  scheduled_for_infection <- union(
    api$get_scheduled(events$infection),
    api$get_scheduled(events$asymptomatic_infection)
  )
  to_infect <- setdiff(
    clinical_infections,
    c(scheduled_for_infection, treated)
  )
  to_infect_asym <- setdiff(
    infections,
    c(scheduled_for_infection, clinical_infections)
  )

  if(length(to_infect) > 0) {
    api$schedule(events$infection, to_infect, parameters$de)
  }

  if(length(to_infect_asym) > 0) {
    api$schedule(events$asymptomatic_infection, to_infect_asym, parameters$de)
  }
}

# =================
# Utility functions
# =================

# Implemented from Winskill 2017 - Supplementary Information page 3
# and Griffin et al 2010 S1 page 7
eir <- function(
  age,
  zeta,
  vector_infectivity,
  parameters
  ) {

  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(zeta, psi)
  .pi * vector_infectivity
}

# Implemented from Griffin et al 2010 S1 page 7
vector_infectivity <- function(infectious_variants, parameters) {
  sum(vnapply(
    seq_along(parameters$blood_meal_rate),
    function(variety) {
      parameters$blood_meal_rate[[variety]] * sum(infectious_variants == variety)
    }
  ))
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

human_pi <- function(zeta, psi) {
 (zeta * psi) / sum(zeta * psi)
}

blood_meal_rate <- function(v, parameters) {
  parameters$blood_meal_rates[v]
}

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

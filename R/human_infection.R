#' @title Simulate malaria infection in humans
#' @description
#' Updates human states and variables to represent asymptomatic/clinical/severe
#' and treated malaria; and resulting boosts in immunity
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param bitten_humans a bitset of bitten humans
#' @param age of each human (timesteps)
#' @param parameters of the model
#' @param timestep current timestep
#' @noRd
simulate_infection <- function(
  variables,
  events,
  bitten_humans,
  age,
  parameters,
  timestep,
  renderer
  ) {
  if (bitten_humans$size() > 0) {
    boost_immunity(
      variables$ib,
      bitten_humans,
      variables$last_boosted_ib,
      timestep,
      parameters$ub
    )
  }

  # Calculate Infected
  infected_humans <- calculate_infections(
    variables,
    bitten_humans,
    parameters,
    renderer,
    timestep
  )

  if (infected_humans$size() > 0) {
    boost_immunity(
      variables$ica,
      infected_humans,
      variables$last_boosted_ica,
      timestep,
      parameters$uc
    )
    boost_immunity(
      variables$id,
      infected_humans,
      variables$last_boosted_id,
      timestep,
      parameters$ud
    )
  }

  clinical_infections <- calculate_clinical_infections(
    variables,
    infected_humans,
    parameters,
    renderer,
    timestep
  )

  update_severe_disease(
    timestep,
    infected_humans,
    variables,
    parameters,
    renderer
  )

  treated <- calculate_treated(
    variables,
    clinical_infections,
    parameters,
    timestep,
    renderer
  )

  renderer$render('n_infections', infected_humans$size(), timestep)

  schedule_infections(
    variables,
    clinical_infections,
    treated,
    infected_humans,
    parameters,
    timestep
  )
}

#' @title Calculate overall infections for bitten humans
#' @description
#' Sample infected humans given prophylaxis and vaccination
#' @param variables a list of all of the model variables
#' @param bitten_humans bitset of bitten humans
#' @param parameters model parameters
#' @param renderer model render object
#' @param timestep current timestep
#' @noRd
calculate_infections <- function(
  variables,
  bitten_humans,
  parameters,
  renderer,
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
    prophylaxis[medicated] <- weibull_survival(
      timestep - drug_time,
      parameters$drug_prophylaxis_shape[drug],
      parameters$drug_prophylaxis_scale[drug]
    )
  }

  # calculate vaccine efficacy
  vaccine_efficacy <- rep(0, length(source_vector))
  vaccine_times <- variables$pev_timestep$get_values(source_vector)
  vaccinated <- vaccine_times > -1
  pev_profile <- variables$pev_profile$get_values(source_vector)
  pev_profile <- pev_profile[vaccinated]
  if (length(vaccinated) > 0) {
    antibodies <- calculate_pev_antibodies(
      timestep - vaccine_times[vaccinated],
      exp(sample_pev_param(pev_profile, parameters$pev_profiles, 'cs')),
      invlogit(sample_pev_param(pev_profile, parameters$pev_profiles, 'rho')),
      exp(sample_pev_param(pev_profile, parameters$pev_profiles, 'ds')),
      exp(sample_pev_param(pev_profile, parameters$pev_profiles, 'dl')),
      parameters
    )
    vmax <- vnapply(parameters$pev_profiles, function(p) p$vmax)
    beta <- vnapply(parameters$pev_profiles, function(p) p$beta)
    alpha <- vnapply(parameters$pev_profiles, function(p) p$alpha)
    vaccine_efficacy[vaccinated] <- calculate_pev_efficacy(
      antibodies,
      vmax[pev_profile],
      beta[pev_profile],
      alpha[pev_profile]
    )
  }

  prob <- b * (1 - prophylaxis) * (1 - vaccine_efficacy)
  infected <- bitset_at(source_humans, bernoulli_multi_p(prob))

  incidence_renderer(
    variables$birth,
    renderer,
    infected,
    source_humans,
    prob,
    'inc_',
    parameters$incidence_rendering_min_ages,
    parameters$incidence_rendering_max_ages,
    timestep
  )

  infected
}

#' @title Calculate clinical infections
#' @description
#' Sample clinical infections from all infections
#' @param variables a list of all of the model variables
#' @param infections bitset of infected humans
#' @param parameters model parameters
#' @param renderer model render
#' @param timestep current timestep
#' @noRd
calculate_clinical_infections <- function(
  variables,
  infections,
  parameters,
  renderer,
  timestep
  ) {
  ica <- variables$ica$get_values(infections)
  icm <- variables$icm$get_values(infections)
  phi <- clinical_immunity(ica, icm, parameters)
  clinical_infections <- bitset_at(infections, bernoulli_multi_p(phi))
  incidence_renderer(
    variables$birth,
    renderer,
    clinical_infections,
    infections,
    phi,
    'inc_clinical_',
    parameters$clinical_incidence_rendering_min_ages,
    parameters$clinical_incidence_rendering_max_ages,
    timestep
  )
  clinical_infections
}

#' @title Calculate severe infections
#' @description
#' Sample severely infected humans from clinically infected
#' @param timestep current timestep
#' @param infections indices of all infected humans
#' @param variables a list of all of the model variables
#' @param parameters model parameters
#' @param renderer model outputs
#' @noRd
update_severe_disease <- function(
  timestep,
  infections,
  variables,
  parameters,
  renderer
  ) {
  age <- get_age(variables$birth$get_values(infections), timestep)
  iva <- variables$iva$get_values(infections)
  ivm <- variables$ivm$get_values(infections)
  theta <- severe_immunity(
    age,
    iva,
    ivm,
    parameters
  )
  develop_severe <- bernoulli_multi_p(theta)
  severe_infections <- bitset_at(infections, develop_severe)
  incidence_renderer(
    variables$birth,
    renderer,
    severe_infections,
    infections,
    theta,
    'inc_severe_',
    parameters$severe_incidence_rendering_min_ages,
    parameters$severe_incidence_rendering_max_ages,
    timestep
  )
  boost_immunity(
    variables$iva,
    infections,
    variables$last_boosted_iva,
    timestep,
    parameters$uv
  )
}

#' @title Calculate treated humans
#' @description
#' Sample treated humans from the clinically infected
#' @param variables a list of all of the model variables
#' @param clinical_infections a bitset of clinically infected humans
#' @param parameters model parameters
#' @param timestep the current timestep
#' @param renderer simulation renderer
#' @noRd
calculate_treated <- function(
    variables,
    clinical_infections,
    parameters,
    timestep,
    renderer
) {
  
  # Render the number of individuals clinically infected in the current timestep
  renderer$render('n_clin_infected', clinical_infections$size(), timestep)
  
  # If the number of clinical infections is 0, return an empty Bitset of treated individuals:
  if(clinical_infections$size() == 0) {
    return(individual::Bitset$new(parameters$human_population))
  }
  
  # Get the treatment coverages and sum them to get total coverage in the current time step:
  treatment_coverages <- get_treatment_coverages(parameters, timestep)
  ft <- sum(treatment_coverages)
  
  # If coverage is no-existent, return a blank Bitset of treated individuals:
  if (ft == 0) {
    return(individual::Bitset$new(parameters$human_population))
  }
  
  # Render the total coverage in the current time step:
  renderer$render('ft', ft, timestep)
  
  # Get the clinically infected individuals who receive treatment and the number of them:
  seek_treatment <- sample_bitset(clinical_infections, ft)
  n_treat <- seek_treatment$size()
  
  # Render the number of people receiving treatment:
  renderer$render('n_treated', n_treat, timestep)
  
  # Assign each individual a drug:
  drugs <- as.numeric(parameters$clinical_treatment_drugs[
    sample.int(
      length(parameters$clinical_treatment_drugs),
      n_treat,
      prob = treatment_coverages,
      replace = TRUE
    )
  ])
  
  # If antimalarial resistance is simulated:
  if(parameters$antimalarial_resistance == TRUE) {
    
    # Render the number of people with susceptible infection residence times at the start of the time step:
    renderer$render('dt_susceptible', sum(variables$dt$get_values() == parameters$dt), timestep)
    # Render the number of people with slow parasite clearance residence times at the start of the time step:
    renderer$render('dt_spc', sum(variables$dt$get_values() == parameters$dt_slow_parasite_clearance), timestep)
    
    # Get each individuals resistance parameters:
    resistance_parameters <- get_antimalarial_resistance_parameters(
      parameters = parameters,
      drugs = drugs,
      timestep = timestep
    )
    
    # Calculate each individuals probability of experiencing slow parasite clearance:
    slow_parasite_clearance_prob <- resistance_parameters$artemisinin_resistance_proportion * resistance_parameters$slow_parasite_clearance_probability  
    
    # Determine which individuals get slow parasite clearance:
    slow_parasite_clearance_indices <- bernoulli_multi_p(slow_parasite_clearance_prob)
    
    # Get Bitset of individuals who will experience slow parasite clearance:
    slow_parasite_clearance_individuals <- bitset_at(seek_treatment, slow_parasite_clearance_indices)
    
    # Create a Bitset of the non-SPC, still seeking treatment people:
    susceptible_to_treatment <- seek_treatment$copy()$set_difference(slow_parasite_clearance_individuals)
    
    # Get indices of individuals who are successfully treated:
    successful <- bernoulli_multi_p(parameters$drug_efficacy[drugs])
    
    # Get the indices of the successfully treated individuals:
    treated_index <- bitset_at(seek_treatment, successful)
    
    # Render the number of people who experience slow parasite clearance:
    renderer$render('n_slow_parasite_clearance', slow_parasite_clearance_individuals$size(), timestep)
    
    # Render the number of people who do not experience slow parasite clearance:
    renderer$render('n_susceptible_SPC', susceptible_to_treatment$size(), timestep)
    
  } else {
    
    # In absence of resistance, calculate the number of successfully treated people as normal:
    successful <- bernoulli_multi_p(parameters$drug_efficacy[drugs])
    # Retrieve the indices of the successfully treated individuals:
    treated_index <- bitset_at(seek_treatment, successful)
    
  }
  
  # Render the number of successfully treated individuals in the current time step:
  renderer$render('n_successfully_treated', treated_index$size(), timestep)
  
  # Update those who have been treated
  if (treated_index$size() > 0) {
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
    if(parameters$antimalarial_resistance) {
      variables$dt$queue_update(
        parameters$dt,
        susceptible_to_treatment
      )
      variables$dt$queue_update(
        resistance_parameters$dt_slow_parasite_clearance[slow_parasite_clearance_indices],
        slow_parasite_clearance_individuals
      )
    }
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
#' @param parameters model parameters
#' @noRd
schedule_infections <- function(
  variables,
  clinical_infections,
  treated,
  infections,
  parameters,
  timestep
  ) {
  included <- treated$not(TRUE)

  to_infect <- clinical_infections$and(included)
  to_infect_asym <- clinical_infections$copy()$not(TRUE)$and(infections)$and(
    included
  )

  if(to_infect$size() > 0) {
    update_infection(
      variables$state,
      'D',
      variables$infectivity,
      parameters$cd,
      to_infect
    )
  }

  if(to_infect_asym$size() > 0) {
    update_to_asymptomatic_infection(
      variables,
      parameters,
      timestep,
      to_infect_asym
    )
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
  acquired_immunity[acquired_immunity > 0] <- acquired_immunity[acquired_immunity > 0] + 0.5
  
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
  acquired_immunity[acquired_immunity > 0] <- acquired_immunity[acquired_immunity > 0] + 0.5
  
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

# Implemented from Winskill 2017 - Supplementary Information page 4
blood_immunity <- function(ib, parameters) {
  ib[ib > 0] <- ib[ib > 0] + 0.5
  
  parameters$b0 * (
    parameters$b1 +
      (1 - parameters$b1) /
      (1 + (ib / parameters$ib0) ** parameters$kb)
  )
}

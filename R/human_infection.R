#' @title Simulate malaria infection in humans
#' @description
#' This function ends with the assignment of rates of infection to the competing
#' hazard resolution object and boosts immunity given infectious bites.
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param bitten_humans a bitset of bitten humans
#' @param age of each human (timesteps)
#' @param parameters of the model
#' @param timestep current timestep
#' @param renderer the model renderer object
#' @param infection_outcome competing hazards object for infection rates
#' @noRd
simulate_infection <- function(
    variables,
    events,
    bitten_humans,
    age,
    parameters,
    timestep,
    renderer,
    infection_outcome
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
  calculate_infections(
    variables,
    bitten_humans,
    parameters,
    renderer,
    timestep,
    infection_outcome
  )
}

#' @title Calculate overall infections for bitten humans
#' @description Infection rates are stored in the infection outcome competing hazards object
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
    timestep,
    infection_outcome
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
  vaccine_times <- variables$last_eff_pev_timestep$get_values(source_vector)
  pev_profile <- variables$pev_profile$get_values(source_vector)
  # get vector of individuals who have received their 3rd dose
  vaccinated <- vaccine_times > -1
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

  ## probability of incidence must be rendered at each timestep
  incidence_probability_renderer(
    variables$birth,
    renderer,
    source_humans,
    prob,
    "inc_",
    parameters$incidence_min_ages,
    parameters$incidence_max_ages,
    timestep
  )
  
  ## capture infection rates to resolve in competing hazards
  infection_outcome$set_rates(
    source_humans,
    prob_to_rate(prob))
}

#' @title Assigns infections to appropriate human states
#' @description
#' Updates human states and variables to represent asymptomatic/clinical/severe
#' and treated malaria; and resulting boosts in immunity
#' @param timestep current timestep
#' @param infected_humans bitset of infected humans
#' @param variables a list of all of the model variables
#' @param renderer model render object
#' @param parameters model parameters
#' @param prob vector of population probabilities of infection
#' @noRd
infection_outcome_process <- function(
    timestep,
    infected_humans,
    variables,
    renderer,
    parameters,
    prob){
  
  incidence_renderer(
    variables$birth,
    renderer,
    infected_humans,
    'inc_',
    parameters$incidence_rendering_min_ages,
    parameters$incidence_rendering_max_ages,
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
    'inc_clinical_',
    parameters$clinical_incidence_rendering_min_ages,
    parameters$clinical_incidence_rendering_max_ages,
    timestep
  )
  incidence_probability_renderer(
    variables$birth,
    renderer,
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
    'inc_severe_',
    parameters$severe_incidence_rendering_min_ages,
    parameters$severe_incidence_rendering_max_ages,
    timestep
  )
  incidence_probability_renderer(
    variables$birth,
    renderer,
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
  
  if(clinical_infections$size() == 0) {
    return(individual::Bitset$new(parameters$human_population))
  }
  
  treatment_coverages <- get_treatment_coverages(parameters, timestep)
  ft <- sum(treatment_coverages)
  
  if (ft == 0) {
    return(individual::Bitset$new(parameters$human_population))
  }
  
  renderer$render('ft', ft, timestep)
  seek_treatment <- sample_bitset(clinical_infections, ft)
  n_treat <- seek_treatment$size()
  renderer$render('n_treated', n_treat, timestep)
  
  drugs <- as.numeric(parameters$clinical_treatment_drugs[
    sample.int(
      length(parameters$clinical_treatment_drugs),
      n_treat,
      prob = treatment_coverages,
      replace = TRUE
    )
  ])
  
  successfully_treated <- calculate_successful_treatments(
    parameters,
    seek_treatment,
    drugs,
    timestep,
    renderer,
    ""
  )
  
  if (successfully_treated$successfully_treated$size() > 0) {
    
    if(parameters$antimalarial_resistance) {
      dt_update_vector <- successfully_treated$dt_spc_combined
    } else {
      dt_update_vector <- parameters$dt
    }
    
    update_infection(
      variables$state,
      'Tr',
      variables$infectivity,
      parameters$cd * parameters$drug_rel_c[successfully_treated$drugs],
      variables$progression_rates,
      1/dt_update_vector,
      successfully_treated$successfully_treated
    )
    
    variables$drug$queue_update(
      successfully_treated$drugs,
      successfully_treated$successfully_treated
    )
    variables$drug_time$queue_update(
      timestep,
      successfully_treated$successfully_treated
    )
  }
  
  successfully_treated$successfully_treated
  
}


#' @title Calculate successfully treated humans
#' @description
#' Sample successful treatments based on drug efficacy and antimalarial resistance
#' @param parameters model parameters
#' @param target bitset of treated humans
#' @param drugs drug index
#' @param timestep the current timestep
#' @param renderer simulation renderer
#' @param int_name the intervention name to use for rendering, use "" for frontline treatment
#' @noRd
calculate_successful_treatments <- function(
    parameters,
    target,
    drugs,
    timestep,
    renderer,
    int_name){
  
  #+++ DRUG EFFICACY +++#
  #+++++++++++++++++++++#
  effectively_treated_index <- bernoulli_multi_p(parameters$drug_efficacy[drugs])
  effectively_treated <- bitset_at(target, effectively_treated_index)
  drugs <- drugs[effectively_treated_index]
  n_drug_efficacy_failures <- target$size() - effectively_treated$size()
  renderer$render(paste0('n_', int_name, 'drug_efficacy_failures'), n_drug_efficacy_failures, timestep)
  
  #+++ ANTIMALARIAL RESISTANCE +++#
  #+++++++++++++++++++++++++++++++#
  if(parameters$antimalarial_resistance) {
    resistance_parameters <- get_antimalarial_resistance_parameters(
      parameters = parameters,
      drugs = drugs, 
      timestep = timestep
    )
    
    #+++ EARLY TREATMENT FAILURE +++#
    #+++++++++++++++++++++++++++++++#
    early_treatment_failure_probability <- resistance_parameters$artemisinin_resistance_proportion * resistance_parameters$early_treatment_failure_probability
    successfully_treated_indices <- bernoulli_multi_p(p = 1 - early_treatment_failure_probability)
    successfully_treated <- bitset_at(effectively_treated, successfully_treated_indices)
    n_early_treatment_failure <- effectively_treated$size() - successfully_treated$size()
    renderer$render(paste0('n_', int_name, 'early_treatment_failure'), n_early_treatment_failure, timestep)
    drugs <- drugs[successfully_treated_indices]
    dt_slow_parasite_clearance <- resistance_parameters$dt_slow_parasite_clearance[successfully_treated_indices]
    
    #+++ SLOW PARASITE CLEARANCE +++#
    #+++++++++++++++++++++++++++++++#
    slow_parasite_clearance_probability <- resistance_parameters$artemisinin_resistance_proportion[successfully_treated_indices] *
      resistance_parameters$slow_parasite_clearance_probability[successfully_treated_indices]
    slow_parasite_clearance_indices <- bernoulli_multi_p(p = slow_parasite_clearance_probability)
    slow_parasite_clearance_individuals <- bitset_at(successfully_treated, slow_parasite_clearance_indices)
    renderer$render(paste0('n_', int_name, 'slow_parasite_clearance'), slow_parasite_clearance_individuals$size(), timestep)
    non_slow_parasite_clearance_individuals <- successfully_treated$copy()$set_difference(slow_parasite_clearance_individuals)
    renderer$render(paste0('n_', int_name, 'successfully_treated'), successfully_treated$size(), timestep)
    dt_slow_parasite_clearance <- dt_slow_parasite_clearance[slow_parasite_clearance_indices]
    
    dt_spc_combined <- rep(parameters$dt, length(successfully_treated_indices))
    dt_spc_combined[slow_parasite_clearance_indices] <- dt_slow_parasite_clearance
    
    successfully_treated_list <- list(
      drugs = drugs,
      successfully_treated = successfully_treated,
      dt_spc_combined = dt_spc_combined)
    
  } else {
    
    successfully_treated <- effectively_treated
    renderer$render(paste0('n_', int_name, 'successfully_treated'), successfully_treated$size(), timestep)
    
    successfully_treated_list <- list(
      drugs = drugs,
      successfully_treated = successfully_treated)
    
  }
  successfully_treated_list
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
      variables$progression_rates,
      1/parameters$dd,
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

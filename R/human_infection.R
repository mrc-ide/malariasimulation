#' @title Simulate malaria infection in humans
#' @description
#' This function ends with the assignment of rates of infection to the competing
#' hazard resolution object.  Boosts immunity given infectious bites.
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
  
  if (bitten_humans$bitten_humans$size() > 0) {
    if(parameters$parasite == "falciparum"){
      boost_immunity(
        variables$ib,
        bitten_humans$bitten_humans,
        variables$last_boosted_ib,
        timestep,
        parameters$ub
      )
    }
    
    # Calculate Infected
    calculate_infection_rates(
      variables,
      bitten_humans,
      parameters,
      renderer,
      timestep,
      infection_outcome
    )
  }
}

#' @title Calculate overall infection rates for bitten humans
#' @description
#' Sample infected humans given prophylaxis and vaccination
#' @param variables a list of all of the model variables
#' @param bitten_humans bitset of bitten humans
#' @param parameters model parameters
#' @param renderer model render object
#' @param timestep current timestep
#' @noRd
calculate_infection_rates <- function(
    variables,
    bitten_humans,
    parameters,
    renderer,
    timestep,
    infection_outcome
) {
  
  if(bitten_humans$bitten_humans$size() == 0){return(bitten_humans$bitten_humans)}
  
  if(parameters$parasite == "falciparum"){
    source_humans <- variables$state$get_index_of(c('S','A','U'))$and(bitten_humans$bitten_humans)
    
    ## P. falciparum models blood immunity
    b <- blood_immunity(variables$ib$get_values(source_humans), parameters)
    
  } else if (parameters$parasite == "vivax"){
    ## Source_humans must include individuals with hypnozoites which may be impacted by prophylaxis/vaccination
    source_humans <- bitten_humans$bitten_humans$copy()$or(variables$hypnozoites$get_index_of(0)$not(TRUE))
    bitten_vector <- bitten_humans$bitten_humans$to_vector()
    ## P. vivax does not model blood immunity
    ## But it must take into account multiple bites per person
    b <- 1-(1-parameters$b)^bitten_humans$n_bites_each
  }
  
  # calculate prophylaxis
  source_vector <- source_humans$to_vector()
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
  
  
  if(parameters$parasite == "falciparum"){
    prob <- b * (1 - prophylaxis) * (1 - vaccine_efficacy)
    infection_rates <- numeric(length = parameters$human_population)
    infection_rates[source_vector] <- prob_to_rate(prob)
    
  } else if (parameters$parasite == "vivax"){
    
    ## Calculated rate of infection for all bitten or with hypnozoites
    infection_rates <- relapse_rates <- variables$hypnozoites$get_values() * parameters$f
    infection_rates[bitten_vector] <- infection_rates[bitten_vector] + prob_to_rate(b)
    relative_rate <- relapse_rates/infection_rates
    relative_rate[is.nan(relative_rate)] <- 0
    
    infection_outcome$set_relative_rates(relative_rate)
    
    ## Get relative rates to get probability bitten over relapse
    prob <- rate_to_prob(infection_rates)
    prob[source_vector] <- prob[source_vector] * (1 - prophylaxis) * (1 - vaccine_efficacy)
    infection_rates <- prob_to_rate(prob)
  }
  
  ## Capture infection rates to resolve in competing hazards
  infection_outcome$set_rates(infection_rates)
  infection_outcome$stash_source_humans(source_humans)
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
infection_process_resolved_hazard <- function(
    timestep,
    infected_humans,
    source_humans,
    variables,
    renderer,
    parameters,
    prob,
    relative_rate = NULL){
  
  renderer$render('n_infections', infected_humans$size(), timestep)
  incidence_renderer(
    variables$birth,
    renderer,
    infected_humans,
    source_humans,
    prob[source_humans$to_vector()],
    'inc_',
    parameters$incidence_rendering_min_ages,
    parameters$incidence_rendering_max_ages,
    timestep
  )
  
  ## For P.f we should only boost SAU infections
  ## For P.v SAUDTr infections all get boosted
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
  
  if(parameters$parasite == "falciparum"){
    
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
    
    patent_infections <- NULL
    
  } else if (parameters$parasite == "vivax"){
    
    # identify bite infections from relapse infections
    relapse_infections <- bitset_at(
      infected_humans,
      bernoulli_multi_p(relative_rate[infected_humans$to_vector()])
    )
    
    renderer$render('n_relapses', relapse_infections$size(), timestep)
    # Render all infections by age
    incidence_renderer(
      variables$birth,
      renderer,
      relapse_infections,
      source_humans,
      prob[source_humans$to_vector()],
      'inc_relapse_',
      parameters$incidence_relapse_rendering_min_ages,
      parameters$incidence_relapse_rendering_max_ages,
      timestep
    )
    
    bite_infections <- infected_humans$copy()$and(relapse_infections$not(inplace = F))
    
    ## Drug prophylaxis may limit formation of new hypnozoite batches
    ls_prophylaxis <- rep(0, bite_infections$size())
    if(length(parameters$drug_hypnozoite_efficacy)>0){
      
      ls_drug <- variables$ls_drug$get_values(bite_infections)
      ls_medicated <- (ls_drug > 0)
      ls_medicated[ls_drug > 0] <- !is.na(parameters$drug_hypnozoite_efficacy[ls_drug])
      
      if (any(ls_medicated)) {
        ls_drug <- ls_drug[ls_medicated]
        ls_drug_time <- variables$ls_drug_time$get_values(bite_infections)[ls_medicated]
        ls_prophylaxis[ls_medicated] <- weibull_survival(
          timestep - ls_drug_time,
          parameters$drug_hypnozoite_prophylaxis_shape[ls_drug],
          parameters$drug_hypnozoite_prophylaxis_scale[ls_drug]
        )
      }
    }
    
    ## All bitten humans with an infectious bite (incorporating prophylaxis) get a new batch of hypnozoites
    if(bite_infections$size()>0){
      new_hypnozoite_batch_formed <- bitset_at(bite_infections, bernoulli_multi_p(1-ls_prophylaxis))
      
      # Cap batches
      new_batch_number <- ifelse(variables$hypnozoites$get_values(new_hypnozoite_batch_formed) == parameters$kmax,
                                 variables$hypnozoites$get_values(new_hypnozoite_batch_formed),
                                 variables$hypnozoites$get_values(new_hypnozoite_batch_formed) + 1)
      
      variables$hypnozoites$queue_update(
        new_batch_number,
        new_hypnozoite_batch_formed$and(bite_infections)
      )
    }
    
    ## Only S and U infections need to be split using the patent infection function
    patent_infections <- calculate_patent_infections(
      variables,
      variables$state$get_index_of(c("S","U"))$and(infected_humans),
      parameters,
      renderer,
      timestep
    )
    
    # Patent level infected S and U, and all A infections to get clinical infections
    clinical_infections <- calculate_clinical_infections(
      variables,
      variables$state$get_index_of("A")$and(infected_humans)$or(patent_infections),
      parameters,
      renderer,
      timestep
    )
  }
  
  treated <- calculate_treated(
    variables,
    clinical_infections,
    parameters,
    timestep,
    renderer
  )
  
  schedule_infections(
    variables,
    patent_infections,
    clinical_infections,
    treated,
    infected_humans,
    parameters,
    timestep
  )
}


#' @title Calculate patent infections (vivax only)
#' @description
#' Sample patent infections from all infections
#' @param variables a list of all of the model variables
#' @param infections bitset of infected humans
#' @param parameters model parameters
#' @param renderer model render
#' @param timestep current timestep
#' @noRd
calculate_patent_infections <- function(
    variables,
    infections,
    parameters,
    renderer,
    timestep
) {
  
  id <- variables$id$get_values(infections)
  idm <- variables$idm$get_values(infections)
  
  philm <- anti_parasite_immunity(
    min = parameters$philm_min, max = parameters$philm_max, a50 = parameters$alm50,
    k = parameters$klm, id = id, idm = idm)
  patent_infections <- bitset_at(infections, bernoulli_multi_p(philm))
  
  ## This incidence will only be the new patent infections...
  # No, I actually want patent infections to be all the new patent infections AND clinical infections.
  # So it would be a subset of total infections.
  incidence_renderer(
    variables$birth,
    renderer,
    patent_infections,
    infections,
    philm,
    'inc_patent_',
    parameters$patent_incidence_rendering_min_ages,
    parameters$patent_incidence_rendering_max_ages,
    timestep
  )
  patent_infections
}

#' @title Calculate clinical infections
#' @description
#' Sample clinical infections from all infections or patent infections (vivax)
#' @param variables a list of all of the model variables
#' @param infections bitset of infected/patent humans
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
  
  # Update liver stage drug effects
  if(length(parameters$drug_hypnozoite_efficacy)>0){
    
    ## When ls stage drug efficacy is NA, this does not results in a success
    hyp_successful <- bernoulli_multi_p(parameters$drug_hypnozoite_efficacy[drugs])
    hyp_treated_index <- bitset_at(seek_treatment, hyp_successful)
    
    if (hyp_treated_index$size() > 0) {
      # Successfully treated hypnozoite batches are removed
      variables$hypnozoites$queue_update(0, hyp_treated_index)
      variables$ls_drug$queue_update(
        drugs[hyp_successful],
        hyp_treated_index
      )
      variables$ls_drug_time$queue_update(
        timestep,
        hyp_treated_index
      )
    }
  }
  
  #+++ DRUG EFFICACY +++#
  #+++++++++++++++++++++#
  effectively_treated_index <- bernoulli_multi_p(parameters$drug_efficacy[drugs])
  effectively_treated <- bitset_at(seek_treatment, effectively_treated_index)
  drugs <- drugs[effectively_treated_index]
  n_drug_efficacy_failures <- n_treat - effectively_treated$size()
  renderer$render('n_drug_efficacy_failures', n_drug_efficacy_failures, timestep)
  
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
    renderer$render('n_early_treatment_failure', n_early_treatment_failure, timestep)
    drugs <- drugs[successfully_treated_indices]
    dt_slow_parasite_clearance <- resistance_parameters$dt_slow_parasite_clearance[successfully_treated_indices]
    
    #+++ SLOW PARASITE CLEARANCE +++#
    #+++++++++++++++++++++++++++++++#
    slow_parasite_clearance_probability <- resistance_parameters$artemisinin_resistance_proportion[successfully_treated_indices] *
      resistance_parameters$slow_parasite_clearance_probability[successfully_treated_indices]
    slow_parasite_clearance_indices <- bernoulli_multi_p(p = slow_parasite_clearance_probability)
    slow_parasite_clearance_individuals <- bitset_at(successfully_treated, slow_parasite_clearance_indices)
    renderer$render('n_slow_parasite_clearance', slow_parasite_clearance_individuals$size(), timestep)
    non_slow_parasite_clearance_individuals <- successfully_treated$copy()$set_difference(slow_parasite_clearance_individuals)
    renderer$render('n_successfully_treated', successfully_treated$size(), timestep)
    dt_slow_parasite_clearance <- dt_slow_parasite_clearance[slow_parasite_clearance_indices]
    
  } else {
    
    successfully_treated <- effectively_treated
    renderer$render('n_successfully_treated', successfully_treated$size(), timestep)
    
  }
  
  if (successfully_treated$size() > 0) {
    variables$state$queue_update("Tr", successfully_treated)
    variables$infectivity$queue_update(
      parameters$cd * parameters$drug_rel_c[drugs],
      successfully_treated
    )
    variables$drug$queue_update(
      drugs,
      successfully_treated
    )
    variables$drug_time$queue_update(
      timestep,
      successfully_treated
    )
    if(parameters$antimalarial_resistance) {
      variables$dt$queue_update(
        parameters$dt,
        non_slow_parasite_clearance_individuals
      )
      variables$dt$queue_update(
        dt_slow_parasite_clearance,
        slow_parasite_clearance_individuals
      )
    }
  }
  
  successfully_treated
  
}

#' @title Schedule infections
#' @description
#' Schedule infections in humans after the incubation period
#' @param events a list of all of the model events
#' @param patent_infections bitset of patent-level infected humans (P.v only)
#' @param clinical_infections bitset of clinically infected humans
#' @param treated bitset of treated humans
#' @param infections bitset of infected humans
#' @param parameters model parameters
#' @noRd
schedule_infections <- function(
    variables,
    patent_infections = NULL,
    clinical_infections,
    treated,
    infections,
    parameters,
    timestep
) {
  
  included <- treated$not(FALSE)
  to_infect <- clinical_infections$and(included)
  
  if(to_infect$size() > 0) {
    update_infection(
      variables$state,
      'D',
      variables$infectivity,
      parameters$cd,
      to_infect
    )
  }
  
  # falciparum infection can result in D or A
  if(parameters$parasite == "falciparum"){
    to_infect_asym <- clinical_infections$copy()$not(FALSE)$and(infections)$and(included)
    
    if(to_infect_asym$size() > 0) {
      # falciparum has age- and immunity-dependent asymptomatic infectivity
      update_to_asymptomatic_infection(
        variables,
        parameters,
        timestep,
        to_infect_asym
      )}
    
    # vivax infection can result in D, A or U
  } else if (parameters$parasite == "vivax"){
    
    ## I think here I can literally just make the appropriate assignments.
    
    to_infect_subpatent <-  variables$state$get_index_of(c('S'))$and(included)$and(infections)$and(patent_infections$not(inplace = F))
    to_infect_asym <-       variables$state$get_index_of(c('S',"U"))$and(included)$and(patent_infections)$and(clinical_infections$not(inplace = F))
    
    if(to_infect_asym$size() > 0) {
      # vivax has constant asymptomatic infectivity
      update_infection(
        variables$state,
        'A',
        variables$infectivity,
        parameters$ca,
        to_infect_asym
      )}
    
    if(to_infect_subpatent$size() > 0) {
      update_infection(
        variables$state,
        'U',
        variables$infectivity,
        parameters$cu,
        to_infect_subpatent
      )
    }
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
    # boost the immunity variable
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
  if(parameters$parasite == "falciparum"){
    # This step only appears in the falciparum model
    acquired_immunity[acquired_immunity > 0] <- acquired_immunity[acquired_immunity > 0] + 0.5
  }
  
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

# Implemented from White et al., 2018 - Supplementary Information
anti_parasite_immunity <- function(min, max, a50, k, id, idm){
  min + (max - min) / (
    1 + ((id + idm) / a50) ** k)
}

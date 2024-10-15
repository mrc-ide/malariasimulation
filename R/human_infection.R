#' @title Simulate malaria infection in humans
#' @description
#' This function ends with the assignment of rates of infection to the competing
#' hazard resolution object and boosts immunity given infectious bites.
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param bitten_humans a bitset of bitten humans
#' @param n_bites_per_person vector of number of bites each person receives (p.v only)
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
    n_bites_per_person,
    age,
    parameters,
    timestep,
    renderer,
    infection_outcome
) {
  
  if (bitten_humans$size() > 0) {
    if(parameters$parasite == "falciparum"){
      
      boost_immunity(
        variables$ib,
        bitten_humans,
        variables$last_boosted_ib,
        timestep,
        parameters$ub
      )
      
      # Calculate Infected
      calculate_falciparum_infections(
        variables,
        bitten_humans,
        parameters,
        renderer,
        timestep,
        infection_outcome
      )
      
    } else if(parameters$parasite == "vivax"){
      
      # Calculate Infected
      calculate_vivax_infections(
        variables,
        bitten_humans,
        n_bites_per_person,
        parameters,
        renderer,
        timestep,
        infection_outcome
      )
    }
  }
  
}

#' @title Calculate overall falciparum infections for bitten humans
#' @description Infection rates are stored in the infection outcome competing hazards object
#' @param variables a list of all of the model variables
#' @param bitten_humans bitset of bitten humans
#' @param parameters model parameters
#' @param renderer model render object
#' @param timestep current timestep
#' @noRd
calculate_falciparum_infections <- function(
    variables,
    bitten_humans,
    parameters,
    renderer,
    timestep,
    infection_outcome
) {
  
  source_humans <- variables$state$get_index_of(c('S','A','U'))$and(bitten_humans)
  
  if(source_humans$size() > 0){
    
    source_vector <- source_humans$to_vector()
    
    # calculate prophylaxis
    prophylaxis <- treatment_prophylaxis_efficacy(
      source_vector,
      variables,
      parameters,
      timestep
    )
    
    # calculate vaccine efficacy
    vaccine_efficacy <- pev_efficacy(
      source_vector,
      variables,
      parameters,
      timestep)
    
    ## p.f models blood immunity
    prob <- blood_immunity(variables$ib$get_values(source_humans), parameters)
    prob <- prob * (1 - prophylaxis) * (1 - vaccine_efficacy)
    infection_rates <- prob_to_rate(prob)
    
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
      infection_rates)
  }
}

#' @title Calculate overall vivax infections for bitten humans
#' @description Infection rates are stored in the infection outcome competing hazards object
#' @param variables a list of all of the model variables
#' @param bitten_humans bitset of bitten humans
#' @param n_bites_per_person vector of number of bites each person receives (p.v only)
#' @param parameters model parameters
#' @param renderer model render object
#' @param timestep current timestep
#' @noRd
calculate_vivax_infections <- function(
    variables,
    bitten_humans,
    n_bites_per_person,
    parameters,
    renderer,
    timestep,
    infection_outcome
) {
  
  ## source_humans must include individuals with hypnozoites which may be impacted by prophylaxis/vaccination
  hypnozoites_humans <- variables$hypnozoites$get_index_of(set = 0)$not(T)
  source_humans <- bitten_humans$copy()$or(hypnozoites_humans)
  
  if(source_humans$size() > 0){
    
    source_vector <- source_humans$to_vector()
    
    # calculate prophylaxis
    prophylaxis <- treatment_prophylaxis_efficacy(
      source_vector,
      variables,
      parameters,
      timestep
    )
    
    # calculate vaccine efficacy
    vaccine_efficacy <- pev_efficacy(
      source_vector,
      variables,
      parameters,
      timestep)
    
    ## p.v does not model blood immunity but must take into account multiple bites per person
    b <- 1 - (1 - parameters$b)^n_bites_per_person[bitten_humans$to_vector()]
    
    ## get infection rates for all bitten or with hypnozoites
    infection_rates <- rep(0, length = source_humans$size())
    infection_rates[bitset_index(source_humans, bitten_humans)] <- prob_to_rate(b)
    
    # Add relapse rates for individuals with hypnozoites
    relative_rates <- NULL
    if(hypnozoites_humans$size()>0){
      relapse_rates <- variables$hypnozoites$get_values(hypnozoites_humans) * parameters$f
      hyp_source_index <- bitset_index(source_humans, hypnozoites_humans)
      infection_rates[hyp_source_index] <- infection_rates[hyp_source_index] + relapse_rates
      # Get and store relative rates for bite/relapse competing hazards resolution
      relative_rates <- relapse_rates/infection_rates[hyp_source_index]
    }
    
    prob <- rate_to_prob(infection_rates) * (1 - prophylaxis) * (1 - vaccine_efficacy)
    infection_rates <- prob_to_rate(prob)
    
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
      infection_rates,
      relative_rates = relative_rates
    )
  }
}

#' @title Calculate protection from infection due to drug prophylaxis
#' @description This function calculates the probability that drug prophylaxis will 
#' protect each individual from infection
#' @param source_vector a vector of individuals with prospective new infections
#' @param variables a list of all of the model variables
#' @param parameters model parameters
#' @param timestep current timestep
#' @noRd
treatment_prophylaxis_efficacy <- function(
    source_vector,
    variables,
    parameters,
    timestep
){
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
  prophylaxis
}

#' @title Calculate protection from infection due to vaccination
#' @description This function calculates the probability that vaccination will 
#' protect each individual from infection
#' @param source_vector a vector of individuals with prospective new infections
#' @param variables a list of all of the model variables
#' @param parameters model parameters
#' @param timestep current timestep
#' @noRd
pev_efficacy <- function(
    source_vector,
    variables,
    parameters,
    timestep
){
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
  vaccine_efficacy
}

#' @title Assigns falciparum infections to appropriate human states
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
falciparum_infection_outcome_process <- function(
    timestep,
    infected_humans,
    variables,
    renderer,
    parameters){
  
  if (infected_humans$size() > 0) {
    
    renderer$render('n_infections', infected_humans$size(), timestep)
    incidence_renderer(
      variables$birth,
      renderer,
      infected_humans,
      'inc_',
      parameters$incidence_rendering_min_ages,
      parameters$incidence_rendering_max_ages,
      timestep
    )
    
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
    
    clinical <- calculate_clinical_infections(
      variables,
      infected_humans,
      parameters,
      renderer,
      timestep
    )
    
    treated <- calculate_treated(
      variables,
      clinical,
      parameters,
      timestep,
      renderer
    )
    
    update_severe_disease(
      timestep,
      infected_humans,
      variables,
      parameters,
      renderer
    )
    
    ## The treated and infected_humans bitsets are re-written so be cautious!
    to_D <- treated$not(FALSE)$and(clinical)
    to_A <- infected_humans$and(clinical$not(FALSE))
    to_U <- NULL
    
    schedule_infections(
      parameters,
      variables,
      timestep,
      to_D,
      to_A,
      to_U
    )
  }
}

#' @title Assigns vivax infections to appropriate human states
#' @description
#' Updates human states and variables to represent asymptomatic/clinical/severe
#' and treated malaria; and resulting boosts in immunity
#' @param timestep current timestep
#' @param infected_humans bitset of infected humans
#' @param variables a list of all of the model variables
#' @param renderer model render object
#' @param parameters model parameters
#' @param relative_rates relative rates of hypnozoite relapse relative to total infection rate
#' @noRd
vivax_infection_outcome_process <- function(
    timestep,
    infected_humans,
    variables,
    renderer,
    parameters,
    relative_rates){
  
  if (infected_humans$size() > 0) {
    
    renderer$render('n_infections', infected_humans$size(), timestep)
    incidence_renderer(
      variables$birth,
      renderer,
      infected_humans,
      'inc_',
      parameters$incidence_rendering_min_ages,
      parameters$incidence_rendering_max_ages,
      timestep
    )
    
    boost_immunity(
      variables$iaa,
      infected_humans,
      variables$last_boosted_iaa,
      timestep,
      parameters$ua
    )
    
    boost_immunity(
      variables$ica,
      infected_humans,
      variables$last_boosted_ica,
      timestep,
      parameters$uc
    )
    
    relapse_bite_infection_hazard_resolution(
      infected_humans,
      relative_rates,
      variables,
      parameters,
      renderer,
      timestep
    )
    
    ## Only S and U infections are considered in generating lm-det infections
    lm_detectable <- calculate_lm_det_infections(
      variables,
      variables$state$get_index_of(c("S","U"))$and(infected_humans),
      parameters
    )
    
    # Lm-detectable level infected S and U, and all A infections may receive clinical infections
    # There is a different calculation to generate clinical infections, based on current infection level
    # LM infections must only pass through the clinical calculation, therefore all "A" infections are included
    # "S" and "U" infections must pass through the lm-detectable calculation prior to and in addition to the clinical
    # calculation. We therefore consider all "A" infections and only the "S" and "U" infections that are now lm-detectable.
    clinical <- calculate_clinical_infections(
      variables,
      variables$state$get_index_of("A")$and(infected_humans)$or(lm_detectable),
      parameters,
      renderer,
      timestep
    )
    
    treated <- calculate_treated(
      variables,
      clinical,
      parameters,
      timestep,
      renderer
    )
    
    ## The infected_humans,lm_detectable and clinical bitsets are re-written so be cautious!
    to_U <- infected_humans$and(lm_detectable$not(F))$and(variables$state$get_index_of(c("S")))
    to_A <- lm_detectable$and(clinical$not(F))
    to_D <- clinical$and(treated$not(F))
    
    schedule_infections(
      parameters,
      variables,
      timestep,
      to_D,
      to_A,
      to_U
    )
  }
}

#' @title Relapse/bite infection competing hazard resolution (p.v only)
#' @description
#' Resolves competing hazards of bite and hypnozoite relapse infections. 
#' For bite infections we increase the batch number and factor in drug prophylaxis.
#' 
#' @param variables a list of all of the model variables
#' @param infected_humans bitset of infected humans
#' @param relative_rates relative rate of relapse infection
#' @param variables model variables
#' @param parameters model parameters
#' @param renderer model renderer
#' @param timestep current timestep
#' @noRd
relapse_bite_infection_hazard_resolution <- function(
    infected_humans,
    relative_rates,
    variables,
    parameters,
    renderer,
    timestep
){
  
  if(variables$hypnozoites$get_index_of(0)$not(T)$and(infected_humans)$size()>0){
    
    hypnozoite_humans <- variables$hypnozoites$get_index_of(0)$not(T)
    potential_relapse_index <- bitset_index(hypnozoite_humans, infected_humans)
    relapse_infections <- bitset_at(
      hypnozoite_humans$and(infected_humans),
      bernoulli_multi_p(relative_rates[potential_relapse_index])
    )
    
    renderer$render('n_relapses', relapse_infections$size(), timestep)
    # render relapse infections by age
    incidence_renderer(
      variables$birth,
      renderer,
      relapse_infections,
      'inc_relapse_',
      parameters$incidence_relapse_rendering_min_ages,
      parameters$incidence_relapse_rendering_max_ages,
      timestep
    )
    
    # get bite infections
    bite_infections <- infected_humans$copy()$and(relapse_infections$not(inplace = F))
    
  } else {
    bite_infections <- infected_humans
  }
  
  ## all bitten humans with an infectious bite (incorporating prophylaxis) get a new batch of hypnozoites
  if(bite_infections$size()>0){
    
    # make sure batches are capped
    current_batches <- variables$hypnozoites$get_values(bite_infections)
    new_batch_number <- pmin(current_batches + 1, parameters$kmax)
    
    variables$hypnozoites$queue_update(
      new_batch_number,
      bite_infections
    )
  }
}

#' @title Calculate light microscopy detectable infections (p.v only)
#' @description
#' Sample light microscopy detectable infections from all infections
#' @param variables a list of all of the model variables
#' @param infections bitset of infected humans
#' @param parameters model parameters
#' @noRd
calculate_lm_det_infections <- function(
    variables,
    infections,
    parameters
) {
  
  iaa <- variables$iaa$get_values(infections)
  iam <- variables$iam$get_values(infections)
  
  philm <- anti_parasite_immunity(
    min = parameters$philm_min, max = parameters$philm_max, a50 = parameters$alm50,
    k = parameters$klm, iaa = iaa, iam = iam)
  
  bitset_at(infections, bernoulli_multi_p(philm))
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
#' @param parameters model parameters
#' @param variables a list of all of the model variables
#' @param timestep current timestep
#' @param to_D bitset of humans to move to state D
#' @param to_A bitset of humans to move to state A
#' @param to_U bitset of humans to move to state U
#' @noRd
schedule_infections <- function(
    parameters,
    variables,
    timestep,
    to_D,
    to_A,
    to_U
) {
  
  if(to_D$size() > 0) {
    update_infection(
      variables$state,
      'D',
      variables$infectivity,
      parameters$cd,
      variables$progression_rates,
      1/parameters$dd,
      to_D
    )
  }
  
  if(to_A$size() > 0) {
    if(parameters$parasite == "falciparum"){
      # p.f has immunity-determined asymptomatic infectivity
      update_to_asymptomatic_infection(
        variables,
        parameters,
        timestep,
        to_A
      )
    } else if (parameters$parasite == "vivax"){
      # p.v has constant asymptomatic infectivity
      update_infection(
        variables$state,
        'A',
        variables$infectivity,
        parameters$ca,
        variables$progression_rates,
        1/parameters$da,
        to_A
      )
    }
  }
  
  if(parameters$parasite == "vivax"){
    # new p.v infections may be pcr-detectable
    if(to_U$size() > 0){
      # p.v pcr-detectable recovery rate is immunity dependent
      update_infection(
        variables$state,
        'U',
        variables$infectivity,
        parameters$cu,
        variables$progression_rates,
        1/anti_parasite_immunity(
          parameters$dpcr_min, parameters$dpcr_max, parameters$apcr50, parameters$kpcr,
          variables$iaa$get_values(to_U),
          variables$iam$get_values(to_U)
        ),
        to_U
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

# Implemented from White et al., 2018 - Supplementary Information
anti_parasite_immunity <- function(min, max, a50, k, iaa, iam){
  min + (max - min) / (
    1 + ((iaa + iam) / a50) ** k)
}
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

  if(parameters$parasite == "falciparum"){

    if (bitten_humans$size() > 0) {
      boost_immunity(
        variables$ib,
        bitten_humans,
        variables$last_boosted_ib,
        timestep,
        parameters$ub
      )}

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

    patent_infections <- calculate_patent_infections(
      variables,
      infected_humans,
      parameters,
      renderer,
      timestep
    )

    clinical_infections <- calculate_clinical_infections(
      variables,
      patent_infections,
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

  if(parameters$parasite == "falciparum"){
    source_humans <- variables$state$get_index_of(c('S','A','U'))$and(bitten_humans)
    source_vector <- source_humans$to_vector()
    ## P. falciparum models blood immunity
    b <- blood_immunity(variables$ib$get_values(source_humans), parameters)

  } else if (parameters$parasite == "vivax"){
    ## All infections must be checked for hypnozoite batch increases
    source_humans <- bitten_humans
    source_vector <- source_humans$to_vector()
    ## P. vivax does not model blood immunity
    b <- parameters$b
  }

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

  if(parameters$parasite == "falciparum"){

    ## All bitten humans with an infectious bite (incorporating prophylaxis and vaccines)
    newly_infected <- bitset_at(source_humans, bernoulli_multi_p(prob))

    # Render new infections caused by bites
    renderer$render('n_infections', newly_infected$size(), timestep)
    incidence_renderer(
      variables$birth,
      renderer,
      newly_infected,
      source_humans,
      prob,
      'inc_',
      parameters$incidence_rendering_min_ages,
      parameters$incidence_rendering_max_ages,
      timestep
    )
  }

  else if(parameters$parasite == "vivax"){

    # Convert prob of infected bite to rate
    rate_bitten <- prob_to_rate(prob)
    rate_bitten_complete <- c(rep(0, parameters$human_population))
    rate_bitten_complete[source_humans$to_vector()] <- rate_bitten

    # Get rates/probs of relapse
    rate_relapse <- variables$hypnozoites$get_values() * parameters$f
    prob_relapse <- rate_to_prob(rate_relapse)

    # Combine and relative
    rate_infection_complete <- rate_bitten_complete + rate_relapse
    relative_rate <- c(rate_bitten_complete/rate_infection_complete)

    # Humans with potential new infections
    with_hypnozoites <- variables$hypnozoites$get_index_of(0)$not(TRUE)
    potential_new_infections_humans <- source_humans$copy()$or(with_hypnozoites)

    # Sum rates of bite infection and relapses
    prob_new_infections <- rate_to_prob(rate_infection_complete)

    # Subset for new infections/bite infections
    newly_infected <- bitset_at(potential_new_infections_humans,
                                bernoulli_multi_p(prob_new_infections[potential_new_infections_humans$to_vector()]))

    newly_bite_infected <- bitset_at(
      newly_infected,
      bernoulli_multi_p(relative_rate[newly_infected$to_vector()]))

    can_increase_batches <- variables$hypnozoites$get_index_of(
      a = 0, b = 9)

    ## Drug prophylaxis may limit formation of new hypnozoite batches
    newly_bite_infected_vector <- newly_bite_infected$to_vector()
    ls_prophylaxis <- rep(0, length(newly_bite_infected_vector))
    if(length(parameters$drug_hypnozoite_efficacy)>0){
      ls_drug <- variables$ls_drug$get_values(newly_bite_infected_vector)
      ls_medicated <- (ls_drug > 0)
      if (any(ls_medicated)) {
        ls_drug <- ls_drug[ls_medicated]
        ls_drug_time <- variables$ls_drug_time$get_values(newly_bite_infected_vector[ls_medicated])
        ls_prophylaxis[ls_medicated] <- weibull_survival(
          timestep - ls_drug_time,
          parameters$drug_hypnozoite_prophylaxis_shape[ls_drug],
          parameters$drug_hypnozoite_prophylaxis_scale[ls_drug]
        )
      }
    }

    ## All bitten humans with an infectious bite (incorporating prophylaxis) get a new batch of hypnozoites
    if(newly_bite_infected$size()>0){
      new_batches_formed <- bitset_at(newly_bite_infected, bernoulli_multi_p(1-ls_prophylaxis))
      variables$hypnozoites$queue_update(
        variables$hypnozoites$get_values(new_batches_formed$and(can_increase_batches)) + 1,
        new_batches_formed$and(newly_bite_infected)
      )
    }

    # Subset to SAU
    SAU_infections <-variables$state$get_index_of(c('S','A','U'))
    newly_infected <- newly_infected$and(SAU_infections)
    newly_bite_infected <- newly_bite_infected$and(SAU_infections)

    # Render new infections caused by bites
    renderer$render('n_new_bite_infections', newly_bite_infected$size(), timestep)
    incidence_renderer(
      variables$birth,
      renderer,
      newly_bite_infected,
      source_humans,
      prob,
      'inc_new_bite_',
      parameters$new_bite_incidence_rendering_min_ages,
      parameters$new_bite_incidence_rendering_max_ages,
      timestep
    )

    # Render relapse infections
    new_relapse_infection <- newly_infected$copy()$and(newly_bite_infected$not(inplace = F))
    renderer$render('n_new_relapse_infections', new_relapse_infection$size(), timestep)
    incidence_renderer(
      variables$birth,
      renderer,
      new_relapse_infection,
      with_hypnozoites,
      prob_relapse,
      'inc_relapse_',
      parameters$relapse_incidence_rendering_min_ages,
      parameters$relapse_incidence_rendering_max_ages,
      timestep
    )
  }

  newly_infected

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

    ## Only ls stage drugs have hypnozoite efficacy, so liver stage drug variable will be consistent in prophylaxis
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

  # Update blood stage drug effects
  successful <- bernoulli_multi_p(parameters$drug_efficacy[drugs])
  treated_index <- bitset_at(seek_treatment, successful)

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
  }
  treated_index
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

    to_infect_subpatent <- infections$and(patent_infections$not(inplace = F))$and(included)$and(variables$state$get_index_of(c('S','U')))
    to_infect_asym <- patent_infections$and(included)$and(clinical_infections$not(inplace = F))

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
anti_parasite_immunity <- function(min, max, a50, k, id, idm){
  min + (max - min) / (
    1 + ((id + idm) / a50) ** k)
}

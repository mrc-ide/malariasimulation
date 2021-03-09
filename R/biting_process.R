#' @title Biting process
#' @description
#' This is the biting process. It results in human and mosquito infection and
#' mosquito death.
#' @param renderer the model renderer object
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param parameters model pararmeters
#' @noRd
create_biting_process <- function(renderer, variables, events, parameters) {
  function(timestep) {
    # Calculate combined EIR
    age <- get_age(variables$birth$get_values(), timestep)

    total_eir <- simulate_bites(
      renderer,
      variables,
      events,
      age,
      parameters,
      timestep
    )

    renderer$render("mean_EIR", mean(total_eir), timestep)

    simulate_infection(
      variables,
      events,
      total_eir,
      age,
      parameters,
      timestep
    )
  }
}

simulate_bites <- function(renderer, variables, events, age, parameters, timestep) {
  total_eir <- 0

  human_infectivity <- variables$infectivity$get_values()
  if (parameters$tbv) {
    human_infectivity <- account_for_tbv(
      timestep,
      human_infectivity,
      variables,
      parameters
    )
  }
  renderer$render('infectivity', mean(human_infectivity), timestep)

  # Calculate pi (the relative biting rate for each human)
  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(variables$zeta$get_values(), psi)
  infectious_index <- variables$mosquito_state$get_index_of('Im')
  susceptible_index <- variables$mosquito_state$get_index_of('Sm')
  adult_index <- variables$mosquito_state$get_index_of('Unborn')$not()

  for (s_i in seq_along(parameters$species)) {
    # Calculate the probabilities of each human being bitten (given
    # interventions)
    species_index <- variables$species$get_index_of(
      parameters$species[[s_i]]
    )
    p_bitten <- prob_bitten(timestep, variables, s_i, parameters)

    Q0 <- parameters$Q0[[s_i]]
    Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
    W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
    renderer$render(paste0('p_repelled_', parameters$species[[s_i]]), Z, timestep)
    renderer$render(paste0('p_feed_survives_', parameters$species[[s_i]]), W, timestep)
    f <- blood_meal_rate(s_i, Z, parameters)

    n_infectious <- infectious_index$copy()$and(species_index)$size()

    lambda <- effective_biting_rate(
      .pi,
      age,
      s_i,
      p_bitten,
      f,
      W,
      parameters
    )

    total_eir <- total_eir + n_infectious * lambda

    susceptible_species_index <- susceptible_index$copy()$and(species_index)

    calculate_mosquito_effects(
      variables,
      human_infectivity,
      lambda,
      events,
      s_i,
      susceptible_species_index,
      species_index$and(adult_index),
      W,
      Z,
      f,
      parameters,
      renderer,
      timestep
    )
  }
  total_eir
}

#' @title Simulate malaria infection in humans
#' @description
#' Updates human states and variables to represent asymptomatic/clinical/severe
#' and treated malaria; and resulting boosts in immunity
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param total_eir a vector of eirs for each human summed across each
#' mosquito species
#' @param age of each human (timesteps)
#' @param parameters of the model
#' @param timestep current timestep
#' @noRd
simulate_infection <- function(
  variables,
  events,
  total_eir,
  age,
  parameters,
  timestep
  ) {
  bitten_humans <- bernoulli_multi_p(total_eir)

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
    timestep
  )

  clinical_infections <- calculate_clinical_infections(
    variables,
    infected_humans,
    parameters,
    timestep
  )

  if (parameters$severe_enabled) {
    update_severe_disease(
      timestep,
      clinical_infections,
      variables,
      infected_humans,
      parameters
    )
  }

  treated <- calculate_treated(
    variables,
    clinical_infections,
    events$recovery,
    parameters,
    timestep
  )

  schedule_infections(
    events,
    clinical_infections,
    treated,
    infected_humans,
    parameters
  )
}

# =================
# Utility functions
# =================

#' @title Calculate the effective biting rate for a species on each human given
#' vector control interventions
#' @description
#' Implemented from Griffin et al 2010 S2 page 6
#' @param .pi relative biting rate for each human
#' @param age of each human (timesteps)
#' @param species to model
#' @param p_bitten the probabilities of feeding given vector controls
#' @param f blood meal rate
#' @param W average probability of a successful bite
#' @param parameters of the model
#' @noRd
effective_biting_rate <- function(.pi, age, species, p_bitten, f, W, parameters) {
  a <- human_blood_meal_rate(f, species, W, parameters)
  a * .pi * p_bitten$prob_bitten / sum(.pi * p_bitten$prob_bitten_survives)
}

human_pi <- function(zeta, psi) {
  (zeta * psi) / sum(zeta * psi)
}

blood_meal_rate <- function(v, z, parameters) {
  gonotrophic_cycle <- get_gonotrophic_cycle(v, parameters)
  interrupted_foraging_time <- parameters$foraging_time / (1 - z)
  1 / (interrupted_foraging_time + gonotrophic_cycle)
}

human_blood_meal_rate <- function(f, v, W, parameters) {
  Q <- 1 - (1 - parameters$Q0[[v]]) / W
  Q * f
}

average_p_repelled <- function(p_repelled, .pi, Q0) {
  Q0 * sum(.pi * p_repelled)
}

average_p_successful <- function(prob_bitten_survives, .pi, Q0) {
  (1 - Q0) + Q0 * sum(.pi *  prob_bitten_survives)
}

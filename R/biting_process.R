#' @title Biting process
#' @description
#' This is the biting process. It results in human and mosquito infection and
#' mosquito death.
#' @param individuals a list of individuals in the model
#' @param states a list of all of the model states
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
create_biting_process <- function(
  individuals,
  states,
  variables,
  events
  ) {
  function(api) {
    parameters <- api$get_parameters()
    timestep <- api$get_timestep()

    # Calculate combined EIR
    age <- get_age(api$get_variable(individuals$human, variables$birth), api$get_timestep())

    total_eir <- simulate_bites(api, individuals, states, variables, age, parameters)
    simulate_infection(api, individuals, states, variables, events, total_eir, age, parameters)
  }
}

simulate_bites <- function(api, individuals, states, variables, age, parameters) {
  total_eir <- 0
  lambda <- rep(NA, length(parameters$blood_meal_rate))

  Sm <- api$get_state(individuals$mosquito, states$Sm)
  Pm <- api$get_state(individuals$mosquito, states$Pm)
  Im <- api$get_state(individuals$mosquito, states$Im)
  species_index <- api$get_variable(
    individuals$mosquito,
    variables$mosquito_variety
  )
  human_infectivity <- api$get_variable(individuals$human, variables$infectivity)

  for (species in seq_along(parameters$blood_meal_rate)) {
    p_bitten <- prob_bitten(
      individuals,
      variables,
      species,
      api,
      parameters
    )

    Z <- mean(p_bitten$prob_repelled)
    f <- blood_meal_rate(species, Z, parameters)

    infectious_species_index <- species_index[Im] == species
    n_infectious <- sum(infectious_species_index)

    species_eir <- eir(
      api,
      individuals$human,
      variables$zeta,
      age,
      species,
      n_infectious,
      p_bitten,
      f,
      parameters
    )

    total_eir <- total_eir + species_eir

    susceptible_species <- Sm[species_index[Sm] == species]
    calculate_mosquito_effects(
      api,
      human_infectivity,
      species_eir,
      individuals,
      states,
      species,
      susceptible_species,
      c(
        susceptible_species,
        Pm[species_index[Pm] == species],
        Im[infectious_species_index]
      ),
      mean(p_bitten$prob_bitten_survives),
      Z,
      f,
      parameters
    )
  }
  total_eir
}

simulate_infection <- function(api, individuals, states, variables, events, total_eir, age, parameters) {
  bitten_humans <- which(bernoulli_multi_p(length(total_eir), total_eir))
  api$render("mean_EIR", mean(total_eir))

  ib <- api$get_variable(individuals$human, variables$ib)
  if (length(bitten_humans) > 0) {
    boost_immunity(
      api,
      individuals$human,
      variables$ib,
      bitten_humans,
      ib[bitten_humans],
      variables$last_boosted_ib,
      api$get_timestep(),
      parameters$ub
    )
  }

  # Calculate Infected
  infected_humans <- calculate_infections(
    api,
    individuals$human,
    states,
    variables,
    bitten_humans,
    ib
  )

  clinical_infections <- calculate_clinical_infections(
    api,
    individuals$human,
    variables,
    infected_humans
  )

  if (parameters$severe_enabled) {
    update_severe_disease(
      api,
      clinical_infections,
      age[clinical_infections],
      individuals$human,
      variables,
      infected_humans
    )
  }

  treated <- calculate_treated(
    api,
    individuals$human,
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

# =================
# Utility functions
# =================

# Implemented from Griffin et al 2010 S2 page 6
eir <- function(api, human, zeta, age, species, n_infectious, p_bitten, f, parameters) {
  a <- human_blood_meal_rate(f, species, mean(p_bitten$prob_bitten_survives), parameters)
  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(api$get_variable(human, zeta), psi)
  infectious_bites <- n_infectious * a * .pi
  omega <- sum(.pi * p_bitten$prob_bitten_survives)
  infectious_bites * p_bitten$prob_bitten / omega
}

human_pi <- function(zeta, psi) {
 (zeta * psi) / sum(zeta * psi)
}

blood_meal_rate <- function(v, z, parameters) {
  f <- parameters$blood_meal_rates[[v]]
  gonotrophic_cycle <- 1 / f - parameters$foraging_time
  interrupted_foraging_time <- parameters$foraging_time / (1 - z)
  1 / (interrupted_foraging_time + gonotrophic_cycle)
}

human_blood_meal_rate <- function(f, v, w, parameters) {
  Q <- 1 - (1 - parameters$Q0[[v]]) / w
  Q * f
}

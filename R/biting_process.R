#' @title Biting process
#' @description
#' This is the biting process. It results in human and mosquito infection and
#' mosquito death.
#' @param renderer the model renderer object
#' @param solvers mosquito ode solvers
#' @param models mosquito ode models
#' @param variables a list of all of the model variables
#' @param events a list of all of the model events
#' @param parameters model pararmeters
#' @param lagged_infectivity a LaggedValue class with historical sums of infectivity
#' @param lagged_eir a LaggedValue class with historical EIRs
#' @noRd
create_biting_process <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  parameters,
  lagged_infectivity,
  lagged_eir
  ) {
  function(timestep) {
    # Calculate combined EIR
    age <- get_age(variables$birth$get_values(), timestep)

    bitten_humans <- simulate_bites(
      renderer,
      solvers,
      models,
      variables,
      events,
      age,
      parameters,
      timestep,
      lagged_infectivity,
      lagged_eir
    )

    simulate_infection(
      variables,
      events,
      bitten_humans,
      age,
      parameters,
      timestep,
      renderer
    )
  }
}

#' @importFrom stats rpois
simulate_bites <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  age,
  parameters,
  timestep,
  lagged_infectivity,
  lagged_eir
  ) {
  bitten_humans <- individual::Bitset$new(parameters$human_population)

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
  zeta <- variables$zeta$get_values()
  .pi <- human_pi(zeta, psi)

  # Get some indices for later
  if (parameters$individual_mosquitoes) {
    infectious_index <- variables$mosquito_state$get_index_of('Im')
    susceptible_index <- variables$mosquito_state$get_index_of('Sm')
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not()
  }

  EIR <- 0

  for (s_i in seq_along(parameters$species)) {
    solver_states <- solver_get_states(solvers[[s_i]])
    p_bitten <- prob_bitten(timestep, variables, s_i, parameters)
    Q0 <- parameters$Q0[[s_i]]
    W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
    Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
    f <- blood_meal_rate(s_i, Z, parameters)
    a <- .human_blood_meal_rate(f, s_i, W, parameters)
    lambda <- effective_biting_rates(a, .pi, p_bitten)

    if (parameters$individual_mosquitoes) {
      species_index <- variables$species$get_index_of(
        parameters$species[[s_i]]
      )$and(adult_index)
      n_infectious <- calculate_infectious_individual(
        s_i,
        variables,
        infectious_index,
        adult_index,
        species_index,
        parameters
      )
    } else {
      n_infectious <- calculate_infectious_compartmental(solver_states)
    }

    lagged_eir[[s_i]]$save(n_infectious * a, timestep)
    species_eir <- lagged_eir[[s_i]]$get(timestep - parameters$de)
    EIR <- EIR + species_eir
    n_bites <- rpois(1, species_eir * mean(psi))
    if (n_bites > 0) {
      bitten_humans$insert(
        fast_weighted_sample(n_bites, lambda)
      )
    }

    infectivity <- lagged_infectivity$get(timestep - parameters$delay_gam)
    lagged_infectivity$save(sum(human_infectivity * .pi), timestep)
    foim <- calculate_foim(a, infectivity)
    renderer$render(paste0('FOIM_', s_i), foim, timestep)
    mu <- death_rate(f, W, Z, s_i, parameters)
    renderer$render(paste0('mu_', s_i), mu, timestep)

    if (parameters$individual_mosquitoes) {
      # update the ODE with stats for ovoposition calculations
      aquatic_mosquito_model_update(
        models[[s_i]],
        species_index$size(),
        f,
        mu
      )

      # update the individual mosquitoes
      susceptible_species_index <- susceptible_index$copy()$and(species_index)

      biting_effects_individual(
        variables,
        foim,
        events,
        s_i,
        susceptible_species_index,
        species_index,
        mu,
        parameters,
        timestep
      )
    } else {
      adult_mosquito_model_update(
        models[[s_i]],
        mu,
        foim,
        solver_states[[ADULT_ODE_INDICES['Sm']]],
        f
      )
    }
  }

  renderer$render('EIR', EIR, timestep)
  renderer$render('n_bitten', bitten_humans$size(), timestep)
  bitten_humans
}


# =================
# Utility functions
# =================

calculate_eir <- function(species, solvers, variables, parameters, timestep) {
  a <- human_blood_meal_rate(species, variables, parameters, timestep)
  infectious <- calculate_infectious(species, solvers, variables, parameters)
  age <- get_age(variables$birth$get_values(), timestep)
  psi <- unique_biting_rate(age, parameters)
  infectious * a * mean(psi)
}

effective_biting_rates <- function(a, .pi, p_bitten) {
  a * .pi * p_bitten$prob_bitten / sum(.pi * p_bitten$prob_bitten_survives)
}

calculate_infectious <- function(species, solvers, variables, parameters) {
  if (parameters$individual_mosquitoes) {
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not()
    return(
      calculate_infectious_individual(
        species,
        variables,
        variables$mosquito_state$get_index_of('Im'),
        adult_index,
        variables$species$get_index_of(
          parameters$species[[species]]
        )$and(adult_index),
        parameters
      )
    )
  }
  calculate_infectious_compartmental(solver_get_states(solvers[[species]]))
}

calculate_infectious_individual <- function(
  species,
  variables,
  infectious_index,
  adult_index,
  species_index,
  parameters
  ) {
  
  infectious_index$copy()$and(species_index)$size()
}

calculate_infectious_compartmental <- function(solver_states) {
  solver_states[[ADULT_ODE_INDICES['Im']]]
}

intervention_coefficient <- function(p_bitten) {
   p_bitten$prob_bitten / sum(p_bitten$prob_bitten_survives)
}

human_pi <- function(zeta, psi) {
  (zeta * psi) / sum(zeta * psi)
}

blood_meal_rate <- function(v, z, parameters) {
  gonotrophic_cycle <- get_gonotrophic_cycle(v, parameters)
  interrupted_foraging_time <- parameters$foraging_time / (1 - z)
  1 / (interrupted_foraging_time + gonotrophic_cycle)
}

human_blood_meal_rate <- function(species, variables, parameters, timestep) {
  age <- get_age(variables$birth$get_values(), timestep)
  psi <- unique_biting_rate(age, parameters)
  zeta <- variables$zeta$get_values()
  p_bitten <- prob_bitten(timestep, variables, species, parameters)
  .pi <- human_pi(zeta, psi)
  Q0 <- parameters$Q0[[species]]
  W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
  Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
  f <- blood_meal_rate(species, Z, parameters)
  .human_blood_meal_rate(f, species, W, parameters)
}

.human_blood_meal_rate <- function(f, v, W, parameters) {
  Q <- 1 - (1 - parameters$Q0[[v]]) / W
  Q * f
}

average_p_repelled <- function(p_repelled, .pi, Q0) {
  Q0 * sum(.pi * p_repelled)
}

average_p_successful <- function(prob_bitten_survives, .pi, Q0) {
  (1 - Q0) + Q0 * sum(.pi *  prob_bitten_survives)
}

# Unique biting rate (psi) for a human of a given age
unique_biting_rate <- function(age, parameters) {
  1 - parameters$rho * exp(- age / parameters$a0)
}

#' @title Calculate the force of infection towards mosquitoes
#'
#' @param a human blood meal rate
#' @param infectivity_sum a sum of infectivity weighted by relative biting rate
#' @noRd
calculate_foim <- function(a, infectivity_sum) {
  a * infectivity_sum
}

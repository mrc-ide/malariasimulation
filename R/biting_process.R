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
#' @noRd
create_biting_process <- function(
  renderer,
  solvers,
  models,
  variables,
  events,
  parameters
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
      timestep
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
  timestep
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
  .pi <- human_pi(variables$zeta$get_values(), psi)

  # Get some indices for later
  if (parameters$hybrid_mosquitoes) {
    infectious_index <- variables$mosquito_state$get_index_of('Im')
    susceptible_index <- variables$mosquito_state$get_index_of('Sm')
    adult_index <- variables$mosquito_state$get_index_of('NonExistent')$not()
    renderer$render('total_M', adult_index$size(), timestep)
  }

  for (s_i in seq_along(parameters$species)) {
    if (!parameters$hybrid_mosquitoes) {
      solver_states <- solver_get_states(solvers[[s_i]])
    }

    if (parameters$hybrid_mosquitoes) {
      species_index <- variables$species$get_index_of(
        parameters$species[[s_i]]
      )$and(adult_index)
      n_infectious <- infectious_index$copy()$and(species_index)$size()
    } else {
      n_infectious <- solver_states[[ADULT_ODE_INDICES['Im']]]
    }

    p_bitten <- prob_bitten(timestep, variables, s_i, parameters)

    Q0 <- parameters$Q0[[s_i]]
    Z <- average_p_repelled(p_bitten$prob_repelled, .pi, Q0)
    W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, Q0)
    renderer$render(
      paste0('p_repelled_', parameters$species[[s_i]]),
      Z,
      timestep
    )
    renderer$render(
      paste0('p_feed_survives_', parameters$species[[s_i]]),
      W,
      timestep
    )
    f <- blood_meal_rate(s_i, Z, parameters)

    lambda <- effective_biting_rate(
      .pi,
      age,
      s_i,
      p_bitten,
      f,
      W,
      parameters
    )

    renderer$render(
      paste0('lambda_', parameters$species[[s_i]]),
      mean(lambda),
      timestep
    )

    n_bites <- rpois(1, n_infectious * sum(lambda))
    if (n_bites > 0) {
      bitten_humans$insert(
        sample.int(
          parameters$human_population,
          n_bites,
          replace = TRUE,
          prob=lambda
        )
      )
    }

    foim <- calculate_foim(human_infectivity, lambda)
    renderer$render(paste0('FOIM_', s_i), foim, timestep)
    mu <- death_rate(f, W, Z, s_i, parameters)
    renderer$render(paste0('mu_', s_i), mu, timestep)

    if (parameters$hybrid_mosquitoes) {
      # update the ODE with stats for ovoposition calculations
      if (parameters$hybrid_mosquitoes) {
        total_M <- species_index$size()
        mosquito_model_update(models[[s_i]], total_M, f, mu)
        renderer$render(paste0('total_M_', s_i), total_M, timestep)
      }

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

  renderer$render('EIR', bitten_humans$size(), timestep)
  bitten_humans
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

# Unique biting rate (psi) for a human of a given age
unique_biting_rate <- function(age, parameters) {
  1 - parameters$rho * exp(- age / parameters$a0)
}

#' @title Calculate the force of infection towards mosquitoes
#'
#' @param human_infectivity a vector of infectivities for each human
#' @param lambda a vector of biting rates for each human
#' @noRd
calculate_foim <- function(human_infectivity, lambda) {
  sum(human_infectivity * lambda)
}

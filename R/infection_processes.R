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
  events,
  odes=NULL
  ) {
  function(api) {
    human <- individuals$human
    parameters <- api$get_parameters()
    timestep <- api$get_timestep()
    source_humans <- api$get_state(
      human,
      states$S,
      states$U,
      states$A
    )

    # Calculate EIR
    age <- get_age(api$get_variable(human, variables$birth), api$get_timestep())

    if (parameters$vector_ode) {
      infectivity <- vector_infectivity_ode(
        odes,
        parameters
      )
    } else {
      source_mosquitos <- api$get_state(individuals$mosquito, states$Im)
      infectivity <- vector_infectivity_ibm(
        api$get_variable(
          individuals$mosquito,
          variables$mosquito_variety,
          source_mosquitos
        ),
        parameters
      )
    }

    epsilon <- eir(
      age,
      api$get_variable(human, variables$xi),
      source_humans,
      infectivity,
      parameters
    )

    api$render("mean_EIR", mean(epsilon))

    bitten_humans <- source_humans[bernoulli(length(source_humans), epsilon)]

    api$render("num_bitten", length(bitten_humans))
    api$render("average_age", mean(age)/365)

    # Calculate Infected
    ib <- api$get_variable(human, variables$ib, bitten_humans)
    b <- blood_immunity(ib, parameters)

    infected_humans <- bitten_humans[bernoulli(length(bitten_humans), b)]

    ica <- api$get_variable(
      human,
      variables$ica,
      infected_humans
    )

    iva <- api$get_variable(
      human,
      variables$iva,
      infected_humans
    )

    icm <- api$get_variable(human, variables$icm, infected_humans)
    phi <- clinical_immunity(
      ica,
      icm,
      parameters
    )

    develop_clinical <- bernoulli(length(infected_humans), phi)

    if (parameters$severe_enabled) {
      theta <- severe_immunity(
        age[infected_humans],
        iva,
        api$get_variable(human, variables$ivm, infected_humans),
        parameters
      )
      develop_severe <- bernoulli(length(develop_clinical), theta)
    }

    # Exclude humans already scheduled for infection
    scheduled_for_infection <- union(
      api$get_scheduled(events$infection),
      api$get_scheduled(events$asymptomatic_infection)
    )
    to_infect <- setdiff(
      infected_humans[develop_clinical],
      scheduled_for_infection
    )
    to_infect_asym <- setdiff(
      infected_humans[!develop_clinical],
      scheduled_for_infection
    )

    last_infected <- api$get_variable(
      human,
      variables$last_infected,
      infected_humans
    )

    # Updates for those who were bitten
    if (length(bitten_humans) > 0) {
      # Boost immunity
      api$queue_variable_update(
        human,
        variables$ib,
        boost_acquired_immunity(
          ib,
          api$get_variable(
            human,
            variables$last_bitten
          )[bitten_humans],
          timestep,
          parameters$ub
        ),
        bitten_humans
      )
      # record last bitten
      api$queue_variable_update(
        human,
        variables$last_bitten,
        timestep,
        bitten_humans
      )

      # Updates for those who were infected
      if (length(infected_humans) > 0) {
        # Boost immunity
        api$queue_variable_update(
          human,
          variables$ica,
          boost_acquired_immunity(
            ica,
            last_infected,
            timestep,
            parameters$uc
          ),
          infected_humans
        )
        api$queue_variable_update(
          human,
          variables$iva,
          boost_acquired_immunity(
            iva,
            last_infected,
            timestep,
            parameters$uv
          ),
          infected_humans
        )
        api$queue_variable_update(
          human,
          variables$id,
          boost_acquired_immunity(
            api$get_variable(human, variables$id, infected_humans),
            last_infected,
            timestep,
            parameters$ud
          ),
          infected_humans
        )
        # record last infected
        api$queue_variable_update(
          human,
          variables$last_infected,
          timestep,
          infected_humans
        )

        # Schedule infection states
        if(length(to_infect) > 0) {
          api$schedule(events$infection, to_infect, parameters$de)

          if(parameters$severe_enabled) {
            is_severe <- rep(0, length(infected_humans))
            is_severe[develop_severe] <- 1
            api$queue_variable_update(
              human,
              variables$is_severe,
              is_severe,
              infected_humans
            )
          }
        }
        if(length(to_infect_asym) > 0) {
          api$schedule(
            events$asymptomatic_infection,
            to_infect_asym,
            parameters$de
          )
        }
      }
    }
  }
}

#' @title Mosquito infection process
#' @description
#' This is the process of infection for mosquitos. It results in a state
#' transition from Sm to Im for infected mosquitos.
#'
#' NOTE: this process will become obsolete when the model is reformulated to
#' model individual mosquitos biting individual humans.
#' @param mosquito the mosquito individual
#' @param human the human individual
#' @param states a list of all of the model states
#' @param variables a list of all of the model variables
#' @param mosquito_infection the mosquito infection event
create_mosquito_infection_process <- function(
  mosquito,
  human,
  states,
  variables,
  mosquito_infection
  ) {
  function(api) {
    lambda <- mosquito_force_of_infection_from_api(
      human,
      states,
      variables,
      api
    )

    for (species in seq_along(lambda)) {
      api$render(paste0('FOIM_', species), lambda[[species]])
    }

    parameters <- api$get_parameters()
    source_mosquitos <- api$get_state(mosquito, states$Sm)
    species <- api$get_variable(
      mosquito,
      variables$mosquito_variety,
      source_mosquitos
    )

    infected = source_mosquitos[
      bernoulli(length(source_mosquitos), lambda[species])
    ]
    api$queue_state_update(mosquito, states$Pm, infected)
    api$schedule(mosquito_infection, infected, parameters$dem)
  }
}

# =================
# Utility functions
# =================

# Implemented from Winskill 2017 - Supplementary Information page 3
# and Griffin et al 2010 S1 page 7
eir <- function(
  age,
  xi,
  source_humans,
  vector_infectivity,
  parameters
  ) {

  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(xi, psi)

  .pi[source_humans] * vector_infectivity
}

# Implemented from Griffin et al 2010 S1 page 7
vector_infectivity_ibm <- function(
  infectious_variants,
  parameters
  ) {
  sum(vnapply(
    seq_along(parameters$blood_meal_rate),
    function(variety) {
      parameters$blood_meal_rate[[variety]] * sum(infectious_variants == variety)
    }
  ))
}

# Implemented from Griffin et al 2010 S1 page 7
vector_infectivity_ode <- function(
  odes,
  parameters
  ) {
  sum(vnapply(
    odes,
    function(ode) {
      mosquito_model_get_states(ode)[[6]] # infected state
    }
  ) * parameters$blood_meal_rates)
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
    1 + (age / (parameters$ad / 365)) ** parameters$gammad
  )
  q <- parameters$d1 + (1 - parameters$d1) / (
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

human_pi <- function(xi, psi) {
 (xi * psi) / sum(xi * psi)
}

#' @title calculate FOIM from API
#' @param human handle for the human individual
#' @param states list of available states
#' @param variables list of available variables
#' @param api the IBM api
mosquito_force_of_infection_from_api <- function(
  human,
  states,
  variables,
  api
  ) {
  parameters <- api$get_parameters()
  age <- get_age(api$get_variable(human, variables$birth), api$get_timestep())
  xi  <- api$get_variable(human, variables$xi)
  a_subset <- api$get_state(human, states$A)
  d_subset <- api$get_state(human, states$D)
  u_subset <- api$get_state(human, states$U)

  a_infectivity <- asymptomatic_infectivity(
    age[a_subset],
    api$get_variable(human, variables$id, a_subset),
    parameters
  )

  infectivity<- c(
    rep(parameters$cd, length(d_subset)),
    a_infectivity,
    rep(parameters$cu, length(u_subset))
  )

  age_subset <- c(age[d_subset], age[a_subset], age[u_subset])
  xi_subset  <- c(xi[d_subset], xi[a_subset], xi[u_subset])

  mosquito_force_of_infection(
    seq_along(parameters$blood_meal_rates),
    age,
    xi,
    infectivity,
    c(d_subset, a_subset, u_subset),
    parameters
  )
}


#' @title calculate FOIM
#' @description  Implemented from Griffin et al 2010 S1 page 7
#' @param v vector of varieties to calculate for
#' @param age vector for complete human population
#' @param xi het vector for complete human population
#' @param infectious_set the indecies for humans which are infectious
#' @param infectivity the infectivity for `infectious_set`
#' @param parameters model parameters
mosquito_force_of_infection <- function(
  v,
  age,
  xi,
  infectivity,
  infectious_set,
  parameters) {

  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(xi, psi)
  mean_infectivity <- sum(.pi[infectious_set] * infectivity)
  blood_meal_rate(v, parameters) * mean_infectivity
}

blood_meal_rate <- function(v, parameters) {
  parameters$blood_meal_rates[v]
}

boost_acquired_immunity <- function(level, last_boosted, timestep, delay) {
  to_boost <- (timestep - last_boosted) > delay | (last_boosted == -1)
  level[to_boost] <- level[to_boost] + 1
  level
}

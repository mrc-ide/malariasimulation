#' @title Infection process
#' @description
#' This is the process of infection for humans. It results in human future state
#' changes for infected humans and boosts in immunity.
#' @param human, the human individual
#' @param mosquito, the mosquito individual
#' @param states, a list of all of the model states
#' @param variables, a list of all of the model variables
create_infection_process <- function(individuals, states, variables, events) {
  human <- individuals$human
  mosquito <- individuals$mosquito
  function(api) {
    parameters <- api$get_parameters()
    timestep <- api$get_timestep()
    source_humans <- api$get_state(
      human,
      states$S,
      states$U,
      states$A
    )

    # Calculate EIR
    source_mosquitos <- api$get_state(mosquito, states$Im)
    age <- api$get_variable(human, variables$age)

    epsilon <- eir(
      age[source_humans],
      api$get_variable(human, variables$xi, source_humans),
      api$get_variable(
        mosquito,
        variables$mosquito_variety,
        source_mosquitos
      ),
      parameters
    )

    api$render("mean_EIR", mean(epsilon))

    bitten_humans <- source_humans[bernoulli(length(source_humans), epsilon)]

    # Calculate Infected
    ib <- api$get_variable(human, variables$ib, bitten_humans)
    b <- blood_immunity(ib, parameters)

    infected_humans <- bitten_humans[!bernoulli(length(bitten_humans), b)]

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

    phi <- clinical_immunity(
      ica,
      api$get_variable(human, variables$icm, infected_humans),
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
      develop_severe <- bernoulli(length(infected_humans), theta)
      symptomatic <- develop_severe | develop_clinical #NOTE: is this AND or OR?
    } else {
      symptomatic <- develop_clinical
    }

    # Exclude humans already scheduled for infection
    scheduled_for_infection <- union(
      api$get_scheduled(events$infection),
      api$get_scheduled(events$asymptomatic_infection)
    )
    to_infect <- setdiff(
      infected_humans[symptomatic],
      scheduled_for_infection
    )

    to_infect_asym <- setdiff(
      infected_humans[!symptomatic],
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
          api$clear_schedule(
            events$subpatent_infection,
            to_infect
          )
          api$clear_schedule(
            events$recovery,
            to_infect
          )

          if(parameters$severe_enabled && length(develop_severe) > 0) {
            api$queue_variable_update(
              human,
              variables$is_severe,
              1,
              infected_humans[develop_severe]
            )
          }
        }
        if(length(to_infect_asym) > 0) {
          api$schedule(
            events$asymptomatic_infection,
            to_infect_asym,
            parameters$de
          )
          api$clear_schedule(
            events$subpatent_infection,
            to_infect_asym
          )
          api$clear_schedule(
            events$recovery,
            to_infect_asym
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
#' @param mosquito, the mosquito individual
#' @param human, the human individual
#' @param states, a list of all of the model states
#' @param variables, a list of all of the model variables
create_mosquito_infection_process <- function(
  mosquito,
  human,
  states,
  variables,
  mosquito_infection
  ) {
  function(api) {
    parameters <- api$get_parameters()
    source_mosquitos <- api$get_state(mosquito, states$Sm)

    age <- api$get_variable(human, variables$age)
    xi  <- api$get_variable(human, variables$xi)
    a_subset <- api$get_state(human, states$A)
    d_subset <- api$get_state(human, states$D)
    u_subset <- api$get_state(human, states$U)

    a_infectivity <- asymptomatic_infectivity(
      age[a_subset],
      api$get_variable(human, variables$id, a_subset),
      parameters
    )

    age_subset <- c(age[d_subset], age[a_subset], age[u_subset])
    xi_subset  <- c(xi[d_subset], xi[a_subset], xi[u_subset])
    infectivity<- c(
      rep(parameters$cd, length(d_subset)),
      a_infectivity,
      rep(parameters$cu, length(u_subset))
    )

    lambda <- mosquito_force_of_infection(
      api$get_variable(
        mosquito,
        variables$mosquito_variety,
        source_mosquitos
      ),
      age_subset,
      xi_subset,
      infectivity,
      parameters
    )

    infected = source_mosquitos[
      bernoulli(length(source_mosquitos), lambda)
    ]
    api$queue_state_update(mosquito, states$Pm, infected)
    api$schedule(mosquito_infection, infected, parameters$dem)
  }
}

# =================
# Utility functions
# =================

# Implemented from Winskill 2017 - Supplementary Information page 3
eir <- function(
  age,
  xi,
  infectious_variants,
  parameters
  ) {

  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(xi, psi)

  infectious_count <- as.data.frame(table(infectious_variants))
  infectious_blood_meal_rate <- blood_meal_rate(
    infectious_count$infectious_variants,
    parameters
  )

  epsilon0 <- .pi * sum(infectious_blood_meal_rate * infectious_count$Freq)
  epsilon0 * xi * psi
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
# NOTE: I believe there is a typo on equation (9) and there should be a + 1 on
# the denominator
asymptomatic_infectivity <- function(age, immunity, parameters) {
  fd <- 1 - (1 - parameters$fd0) / (
    1 + (age / parameters$ad) ** parameters$gammad
  )
  q <- parameters$d1 + (1 - parameters$dmin) / (
    1 + fd * ((1 + immunity) / parameters$id0) ** parameters$kd
  )
  parameters$cu + (parameters$cd + parameters$cu) * q ** parameters$gamma1
}

# Unique biting rate (psi) for a human of a given age
unique_biting_rate <- function(age, parameters) {
  1 - parameters$rho * exp(- age / parameters$a0)
}

# Relative biting rate (xi) drawn from log normal
relative_biting_rate <- function(n, parameters) {
  rlnorm(n, -parameters$sigma**2/2, parameters$sigma**2)
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

# Implemented from Griffin et al 2010 S1 page 7
mosquito_force_of_infection <- function(v, age, xi, infectivity, parameters) {
  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(xi, psi)
  mean_infectivity <- sum(.pi * infectivity)
  blood_meal_rate(v, parameters) * mean_infectivity
}

blood_meal_rate <- function(v, parameters) {
  rates <- c(parameters$av1, parameters$av2, parameters$av3)
  rates[v]
}

boost_acquired_immunity <- function(level, last_boosted, timestep, delay) {
  to_boost <- (timestep - last_boosted) > delay | (last_boosted == -1)
  level[to_boost] <- level[to_boost] + 1
  level
}

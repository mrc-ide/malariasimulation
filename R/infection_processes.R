#' @title Infection process
#' @description
#' This is the process of infection for humans. It results in human future state
#' changes for infected humans and boosts in immunity.
#' @param human, the human individual
#' @param mosquito, the mosquito individual
#' @param states, a list of all of the model states
#' @param variables, a list of all of the model variables
create_infection_process <- function(human, mosquito, states, variables) {
  function(simulation_frame, timestep, parameters) {
    source_humans <- simulation_frame$get_state(
      human,
      states$S,
      states$U,
      states$A,
      states$D
    )
    next_infection <- simulation_frame$get_variable(
      human,
      variables$infection_schedule
    )
    next_asymptomatic_infection <- simulation_frame$get_variable(
      human,
      variables$asymptomatic_infection_schedule
    )

    # Calculate EIR
    source_mosquitos <- simulation_frame$get_state(mosquito, states$Im)
    age <- simulation_frame$get_variable(human, variables$age)
    ib <- simulation_frame$get_variable(human, variables$ib)

    epsilon <- eir(
      age[source_humans],
      simulation_frame$get_variable(human, variables$xi)[source_humans],
      simulation_frame$get_variable(
        mosquito,
        variables$mosquito_variety
      )[source_mosquitos],
      parameters
    )

    number_of_bites <- round(epsilon)
    bitten_humans <- source_humans[number_of_bites > 0]

    # Calculate Infected
    b <- blood_immunity(ib[bitten_humans], parameters)

    infected_humans <- bitten_humans[
      !bernoulli(
        length(bitten_humans),
        (1 - b) * number_of_bites[number_of_bites > 0]
      )
    ]

    ica <- simulation_frame$get_variable(
      human,
      variables$ica
    )[infected_humans]

    iva <- simulation_frame$get_variable(
      human,
      variables$iva
    )[infected_humans]

    phi <- clinical_immunity(
      ica,
      simulation_frame$get_variable(human, variables$icm)[infected_humans],
      parameters
    )

    theta <- severe_immunity(
      age[infected_humans],
      iva,
      simulation_frame$get_variable(human, variables$ivm)[infected_humans],
      parameters
    )

    develop_clinical <- bernoulli(length(infected_humans), phi)
    develop_severe <- bernoulli(length(infected_humans), theta)
    symptomatic <- develop_severe | develop_clinical

    # Exclude humans already scheduled for infection
    to_infect <- remove_scheduled(
      infected_humans[symptomatic],
      timestep,
      next_infection[symptomatic],
      next_asymptomatic_infection[symptomatic]
    )

    to_infect_asym <- remove_scheduled(
      infected_humans[!symptomatic],
      timestep,
      next_infection[!symptomatic],
      next_asymptomatic_infection[!symptomatic]
    )

    last_infected <- simulation_frame$get_variable(
      human,
      variables$last_infected
    )[infected_humans]

    updates <- list()

    # Updates for those who were bitten
    if (length(bitten_humans) > 0) {
      updates <- c(updates, list(
        # Boost immunity
        individual::VariableUpdate$new(
          human,
          variables$ib,
          boost_acquired_immunity(
            ib[bitten_humans],
            simulation_frame$get_variable(
              human,
              variables$last_bitten
            )[bitten_humans],
            timestep,
            parameters$ub
          ),
          bitten_humans
        ),
        # record last bitten
        individual::VariableUpdate$new(
          human,
          variables$last_bitten,
          timestep,
          bitten_humans
        )
      ))

      # Updates for those who were infected
      if (length(infected_humans) > 0) {
        updates <- c(updates, list(
          # Boost immunity
          individual::VariableUpdate$new(
            human,
            variables$ica,
            boost_acquired_immunity(
              ica,
              last_infected,
              timestep,
              parameters$uc
            ),
            infected_humans
          ),
          individual::VariableUpdate$new(
            human,
            variables$iva,
            boost_acquired_immunity(
              iva,
              last_infected,
              timestep,
              parameters$uv
            ),
            infected_humans
          ),
          individual::VariableUpdate$new(
            human,
            variables$id,
            boost_acquired_immunity(
              simulation_frame$get_variable(human, variables$id)[infected_humans],
              last_infected,
              timestep,
              parameters$ud
            ),
            infected_humans
          ),
          # record last infected
          individual::VariableUpdate$new(
            human,
            variables$last_infected,
            timestep,
            infected_humans
          )
        ))

        # Schedule infection states
        if(length(to_infect) > 0) {
          updates <- c(updates, individual::VariableUpdate$new(
            human,
            variables$infection_schedule,
            timestep + parameters$de,
            to_infect
          ))

          if(length(develop_severe) > 0) {
            updates <- c(updates, individual::VariableUpdate$new(
              human,
              variables$is_severe,
              1,
              infected_humans[develop_severe]
            ))
          }
        }
        if(length(to_infect_asym) > 0) {
          updates <- c(updates, individual::VariableUpdate$new(
            human,
            variables$asymptomatic_infection_schedule,
            timestep + parameters$de,
            to_infect_asym
          ))
        }
      }
    }

    updates
  }
}

create_infection_scheduler <- function(human, A, I, infection_schedule, asymptomatic_infection_schedule) {
  function(simulation_frame, timestep, parameters) {
    infection    <- which(
      simulation_frame$get_variable(human, infection_schedule) == timestep
    )
    asymptomatic <- which(
      simulation_frame$get_variable(human, asymptomatic_infection_schedule) == timestep
    )
    list(
      individual::StateUpdate$new(human, I, infection),
      individual::StateUpdate$new(human, A, asymptomatic)
    )
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
  variables
  ) {
  function(simulation_frame, timestep, parameters) {
    source_mosquitos <- simulation_frame$get_state(mosquito, states$Sm)

    age <- simulation_frame$get_variable(human, variables$age)
    asymptomatic <- simulation_frame$get_state(human, states$A)

    a_infectivity <- asymptomatic_infectivity(
      age[asymptomatic],
      simulation_frame$get_variable(human, variables$id)[asymptomatic],
      parameters
    )

    # Create a dataframe frame with human age, xi and infectivity
    infectivity_frame <- create_infectivity_frame(
      age,
      simulation_frame$get_variable(human, variables$xi),
      list(
        list(simulation_frame$get_state(human, states$D), parameters$cd),
        list(asymptomatic, a_infectivity),
        list(simulation_frame$get_state(human, states$U), parameters$cu)
      )
    )

    lambda <- mosquito_force_of_infection(
      simulation_frame$get_variable(
        mosquito,
        variables$mosquito_variety
      )[source_mosquitos],
      infectivity_frame,
      parameters
    )

    infected = source_mosquitos[
      bernoulli(length(source_mosquitos), lambda)
    ]
    individual::StateUpdate$new(mosquito, states$Im, infected)
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
mosquito_force_of_infection <- function(v, human_frame, parameters) {
  psi <- unique_biting_rate(human_frame$age, parameters)
  .pi <- human_pi(human_frame$xi, psi)
  mean_infectivity <- sum(.pi * human_frame$infectivity)
  blood_meal_rate(v, parameters) * mean_infectivity
}

blood_meal_rate <- function(v, parameters) {
  rates <- c(parameters$av1, parameters$av2, parameters$av3)
  rates[v]
}

create_infectivity_frame <- function(human_age, human_xi, subset_to_param) {
  do.call(
    'rbind',
    lapply(
      subset_to_param,
      function(x) {
        subset <- x[[1]]
        if (length(subset) == 0) {
          return(
            setNames(
              data.frame(
                matrix(ncol = 3, nrow = 0)
              ),
              c('age', 'xi', 'infectivity')
            )
          )
        }
        data.frame(
          age         = human_age[subset],
          xi          = human_xi[subset],
          infectivity = x[[2]]
        )
      }
    )
  )
}

remove_scheduled <- function(subset, timestep, ...) {
  schedules <- list(...)
  for (schedule in schedules) {
    subset <- setdiff(subset, which(schedule >= timestep))
  }
  if (length(subset) == 0) c() else subset
}

boost_acquired_immunity <- function(level, last_boosted, timestep, delay) {
  to_boost <- (timestep - last_boosted) > delay | (last_boosted == -1)
  level[to_boost] <- level[to_boost] + 1
  level
}

# =========
# Processes
# =========

infection_process <- function(simulation_frame, timestep, parameters) {
  source_humans <- simulation_frame$get_state(human, S, U, A, D)
  source_mosquitos <- simulation_frame$get_state(mosquito, Im)
  age_value <- simulation_frame$get_variable(human, age)
  ib_value <- simulation_frame$get_variable(human, ib)[source_humans]

  epsilon <- probability_bitten(
    age_value[source_humans],
    simulation_frame$get_constant(human, xi)[source_humans],
    simulation_frame$get_constant(mosquito, mosquito_variety)[source_mosquitos],
    parameters
  )

  bitten_humans <- source_humans[runif(length(source_humans), 0, 1) > epsilon]

  b <- blood_immunity(ib_value[bitten_humans], parameters)

  bitten_infected_humans <- bitten_humans[runif(length(bitten_humans), 0, 1) > b]

  ica_infected_value <- simulation_frame$get_variable(human, ica)[infected_humans]
  iva_infected_value <- simulation_frame$get_variable(human, iva)[infected_humans]

  phi <- clinical_immunity(
    ica_infected_value,
    simulation_frame$get_variable(human, icm)[infected_humans],
    parameters
  )

  theta <- severe_immunity(
    age_value[infected_humans],
    iva_infected_value,
    simulation_frame$get_variable(human, ivm)[infected_humans],
    parameters
  )

  develop_clinical <- runif(length(infected_humans), 0, 1) > phi
  develop_severe <- runif(length(infected_humans), 0, 1) > theta
  symptomatic <- develop_severe | develop_clinical

  next_infection <- schedule_infection(
    simulation_frame$get_variable(human, infection_schedule)[infected_humans],
    symptomatic,
    timestep + parameters$de
  )

  next_asymptomatic_infection <- schedule_infection(
    simulation_frame$get_variable(
      human,
      asymptomatic_infection_schedule
    )[infected_humans],
    !symptomatic,
    timestep + parameters$de
  )

  last_infected_value <- simulation_frame$get_variable(last_infected)[infected_humans]
  new_ica <- boost_acquired_immunity(
    ica_infected_value,
    last_infected_value,
    timestep,
    parameters$uc
  )

  new_iva <- boost_acquired_immunity(
    iva_infected_value,
    last_infected_value,
    timestep,
    parameters$uv
  )

  new_ib <- boost_acquired_immunity(
    ib_value[bitten_humans],
    simulation_frame$get_variable(last_bitten)[bitten_humans],
    timestep,
    parameters$ub
  )

  list(
    # Boost immunity
    VariableUpdate$new(human, ica, new_ica, infected_humans),
    VariableUpdate$new(human, iva, new_iva, infected_humans),
    VariableUpdate$new(human, ib, new_ib, bitten_humans),
    # Schedule infection states
    VariableUpdate$new(
      human,
      infection_schedule,
      next_infection,
      infected_humans
    ),
    VariableUpdate$new(
      human,
      asymptomatic_infection_schedule,
      next_asymptomatic_infection,
      infected_humans
    ),
    # Record last bitten/infected
    VariableUpdate$new(human, last_bitten, timestep, bitten_humans),
    VariableUpdate$new(human, last_infected, timestep, bitten_infected_humans)
  )
}

scheduled_infections <- function(simulation_frame, timestep, parameters) {
  infection    <- which(
    simulation_frame$get_variable(human, infection_schedule) == timestep
  )
  asymptomatic <- which(
    simulation_frame$get_variable(human, asymptomatic_infection_schedule) == timestep
  )
  list(
    StateUpdate$new(human, I, infection),
    StateUpdate$new(human, A, asymptomatic)
  )
}

mosquito_infection_process <- function(simulation_frame, timestep, parameters) {
  source_mosquitos <- simulation_frame$get_state(mosquito, Sm)

  # Create a dataframe frame with human age, xi and infectivity
  infectivity_frame <- create_infectivity_frame(
    simulation_frame$get_variable(human, age),
    simulation_frame$get_constant(human, xi),
    list(
      list(simulation_frame$get_state(human, D), parameters$cd),
      list(simulation_frame$get_state(human, Treated), parameters$ct),
      list(simulation_frame$get_state(human, A), parameters$ca),
      list(simulation_frame$get_state(human, U), parameters$cu)
    )
  )

  lambda <- mosquito_force_of_infection(
    simulation_frame$get_constant(mosquito, mosquito_variety)[source_mosquitos],
    infectivity_frame,
    parameters
  )
  infected = source_mosquitos[
    runif(length(source_mosquitos), 0, 1) > lambda
  ]
  StateUpdate$new(mosquito, Im, infected)
}

# =================
# Utility functions
# =================

# Implemented from Winskill 2017 - Supplementary Information page 3
probability_bitten <- function(
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
      (acquired_immunity+maternal_immunity) / parameters$iv0) ** parameters$kv
    )
  )
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

schedule_infection <- function(current_schedule, subset, next_event) {
  new_schedule <- current_schedule
  new_schedule[intersect(subset, which(new_schedule == -1))] <- next_event
  new_schedule[intersect(subset, which(new_schedule != -1))] <- pmin(
    new_schedule[subset],
    next_event
  )
  new_schedule
}

boost_acquired_immunity <- function(level, last_boosted, timestep, delay) {
  to_boost <- (timestep - last_boosted) > delay
  level[to_boost] <- level[to_boost] + 1
  level
}

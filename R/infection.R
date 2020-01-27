# =========
# Processes
# =========

infection_process <- function(simulation_frame, timestep, parameters) {
  source_humans <- simulation_frame$get_state(human, S, U, A, D)
  next_infection <- simulation_frame$get_variable(human, infection_schedule)
  next_asymptomatic_infection <- simulation_frame$get_variable(human, infection_schedule)

  # Calculate EIR
  source_mosquitos <- simulation_frame$get_state(mosquito, Im)
  age_value <- simulation_frame$get_variable(human, age)
  ib_value <- simulation_frame$get_variable(human, ib)[source_humans]

  epsilon <- probability_bitten(
    age_value[source_humans],
    simulation_frame$get_constant(human, xi)[source_humans],
    simulation_frame$get_constant(mosquito, mosquito_variety)[source_mosquitos],
    parameters
  )

  bitten_humans <- source_humans[uniform_gt(length(source_humans), epsilon)]

  # Calculate Infected
  b <- blood_immunity(ib_value[bitten_humans], parameters)

  infected_humans <- bitten_humans[uniform_gt(length(bitten_humans), b)]

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

  develop_clinical <- uniform_gt(length(infected_humans), phi)
  develop_severe <- uniform_gt(length(infected_humans), theta)
  symptomatic <- develop_severe | develop_clinical

  # Exclude humans already scheduled for infection
  to_infect <- remove_scheduled(
    infected_humans[symptomatic],
    timestep,
    next_infection,
    next_asymptomatic_infection
  )

  to_infect_asym <- remove_scheduled(
    infected_humans[!symptomatic],
    timestep,
    next_infection,
    next_asymptomatic_infection
  )

  last_infected_value <- simulation_frame$get_variable(
    human,
    last_infected
  )[infected_humans]

  list(
    # Boost immunity
    individual::VariableUpdate$new(
      human,
      ica,
      boost_acquired_immunity(
        ica_infected_value,
        last_infected_value,
        timestep,
        parameters$uc
      ),
      infected_humans
    ),
    individual::VariableUpdate$new(
      human,
      iva,
      boost_acquired_immunity(
        iva_infected_value,
        last_infected_value,
        timestep,
        parameters$uv
      ),
      infected_humans
    ),
    individual::VariableUpdate$new(
      human,
      id,
      boost_acquired_immunity(
        simulation_frame$get_variable(human, id)[infected_humans],
        last_infected_value,
        timestep,
        parameters$ud
      ),
      infected_humans
    ),
    individual::VariableUpdate$new(
      human,
      ib,
      boost_acquired_immunity(
        ib_value[bitten_humans],
        simulation_frame$get_variable(human, last_bitten)[bitten_humans],
        timestep,
        parameters$ub
      ),
      bitten_humans
    ),
    # Schedule infection states
    individual::VariableUpdate$new(
      human,
      infection_schedule,
      timestep + parameters$de,
      to_infect
    ),
    individual::VariableUpdate$new(
      human,
      asymptomatic_infection_schedule,
      timestep + parameters$de,
      to_infect_asym
    ),
    # Record last bitten/infected/is_severe
    individual::VariableUpdate$new(human, last_bitten, timestep, bitten_humans),
    individual::VariableUpdate$new(human, last_infected, timestep, infected_humans),
    individual::VariableUpdate$new(human, is_severe, 1, infected_humans[develop_severe])
  )
}

treatment_recovery_process <- function(simulation_frame, timestep, parameters) {
  source_individuals <- simulation_frame$get_state(human, Treated)
  target_individuals <- source_individuals[
    runif(length(source_individuals), 0, 1) > parameters$rt
  ]
  list(
    individual::StateUpdate$new(human, S, target_individuals),
    individual::VariableUpdate$new(human, is_severe, 0, target_individuals)
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
    individual::StateUpdate$new(human, I, infection),
    individual::StateUpdate$new(human, A, asymptomatic)
  )
}

mosquito_infection_process <- function(simulation_frame, timestep, parameters) {
  source_mosquitos <- simulation_frame$get_state(mosquito, Sm)

  age_value <- simulation_frame$get_variable(human, age)
  asymptomatic <- simulation_frame$get_state(human, A)

  a_infectivity <- asymptomatic_infectivity(
    age_value[asymptomatic],
    simulation_frame$get_variable(human, id)[asymptomatic],
    parameters
  )

  # Create a dataframe frame with human age, xi and infectivity
  infectivity_frame <- create_infectivity_frame(
    age_value,
    simulation_frame$get_constant(human, xi),
    list(
      list(simulation_frame$get_state(human, D), parameters$cd),
      list(simulation_frame$get_state(human, Treated), parameters$ct),
      list(asymptomatic, a_infectivity),
      list(simulation_frame$get_state(human, U), parameters$cu)
    )
  )

  lambda <- mosquito_force_of_infection(
    simulation_frame$get_constant(mosquito, mosquito_variety)[source_mosquitos],
    infectivity_frame,
    parameters
  )

  infected = source_mosquitos[
    uniform_gt(length(source_mosquitos), lambda)
  ]
  individual::StateUpdate$new(mosquito, Im, infected)
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

# =========
# Processes
# =========

infection_process <- function(simulation_frame, timestep, parameters) {
  source_humans <- simulation_frame$get_state(human, S, U, A)
  source_mosquitos <- simulation_frame$get_state(mosquito, Im)

  lambda <- force_of_infection(
    simulation_frame$get_variable(human, age)[source_humans],
    simulation_frame$get_constant(human, xi)[source_humans],
    simulation_frame$get_constant(mosquito, mosquito_variety)[source_mosquitos],
    simulation_frame$get_variable(human, ib)[source_humans],
    parameters
  )

  infected_humans <- source_humans[runif(length(source_humans), 0, 1) > lambda]

  phi <- immunity(
    simulation_frame$get_variable(human, ica)[infected_humans],
    simulation_frame$get_variable(human, icm)[infected_humans],
    parameters
  )

  symptomatic <- runif(length(infected_humans), 0, 1) > phi

  list(
    StateUpdate$new(human, I, infected_humans[symptomatic]),
    StateUpdate$new(human, A, infected_humans[!symptomatic]),
    VariableUpdate$new(
      human,
      last_bitten,
      timestep,
      infected_humans
    ),
    VariableUpdate$new(
      human,
      last_infected,
      timestep,
      infected_humans[symptomatic]
    )
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
# Calculate the unique biting rate (psi) from age
# Calculate the mean EIR (epsilon0) from time
# Sample the relative biting rate (xi) from a normal distribution
# Calculate immunity level (b)
force_of_infection <- function(
  age,
  xi,
  infectious_variants,
  ib,
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
  b <- infection_probability(ib, parameters)
  epsilon0 * xi * b * psi
}

# Implemented from Winskill 2017 - Supplementary Information page 4
# Calculate acquired immunity from last_bitten
# Calculate and maternal immunity from age
# Then calculate immunity using parameters
immunity <- function(acquired_immunity, maternal_immunity, parameters) {
  parameters$theta0 * (
    parameters$theta1 +
      (1 - parameters$theta1) /
      1 + ((acquired_immunity + maternal_immunity) / parameters$ic0)
      ** parameters$kc
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

infection_probability <- function(ib, parameters) {
  parameters$b0 + (parameters$b1 - parameters$b0) /
    (1 + (ib / parameters$ib0)**parameters$kb)
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

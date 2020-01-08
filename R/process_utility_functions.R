# ===============
# Human functions
# ===============

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
    infectious_count$infectious_count,
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
  parameters$phi0 * (
    parameters$phi1 +
      (1 - parameters$phi1) /
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

immunity_decay <- function(level, last_timestep, timestep, rate, delay) {
  boost <- (timestep - last_timestep) == delay
  level[boost] <- level[boost] + 1
  level[!boost] <- level[!boost] - rate * level[!boost]
  level
}

# ==================
# Mosquito functions
# ==================

# Implemented from Griffin et al 2010 S1 page 7
mosquito_force_of_infection <- function(v, human_frame, parameters) {
  psi <- unique_biting_rate(human_frame$age, parameters)
  .pi <- human_pi(human_frame$xi, psi)
  mean_infectivity <- sum(.pi * human_frame$infectivity)
  blood_meal_rate(v, parameters) * mean_infectivity
}

blood_meal_rate <- function(v, parameters) {
  rate <- v
  rate[rate == 1] <- parameters$av1
  rate[rate == 2] <- parameters$av2
  rate[rate == 3] <- parameters$av3
  rate
}

# =================
# General functions
# =================
create_fixed_probability_state_change_process <- function(i, from, to, rate) {
  function (simulation_frame, timestep, parameters) {
    source_individuals <- simulation_frame$get_state(i, from)
    target_individuals <- source_individuals[
      runif(length(source_individuals), 0, 1) > rate
    ]
    StateUpdate$new(i, to, target_individuals)
  }
}

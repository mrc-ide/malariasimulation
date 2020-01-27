#' This is the process for human mortality, it defines which humans die from
#' natural causes and severe infection and replaces dead individuals with
#' newborns.
#'
#' @param simulation_frame, the current state of the simulation
#' @param timestep, the current timestep
#' @param parameters, the model parameters
mortality_process <- function(simulation_frame, timestep, parameters) {
  age_value <- simulation_frame$get_variable(human, age)

  natural_deaths <- died_naturally(
    age_value,
    parameters$mortality_probability_table
  )

  severe_deaths <- died_from_severe(
    which(simulation_frame$get_variable(human, is_severe) == 1),
    simulation_frame$get_state(human, D),
    simulation_frame$get_state(human, Treated),
    parameters$ftv,
    parameters$v
  )

  died <- c(natural_deaths, severe_deaths)

  # workaround index by numeric(0) bug
  if (length(died) == 0) {
    died <- NULL
  }

  # Calculate new maternal immunities
  groups <- simulation_frame$get_variable(human, xi_group)
  sampleable <- age_value >= 15 | age_value <= 35
  icm_values <- simulation_frame$get_variable(human, icm)
  ivm_values <- simulation_frame$get_variable(human, ivm)
  mothers <- sample_mothers(sampleable, died, groups)
  birth_icm <- icm_values[mothers] * parameters$pm
  birth_ivm <- ivm_values[mothers] * parameters$pm

  list(
    individual::VariableUpdate$new(human, age, 0, died),
    individual::VariableUpdate$new(human, last_bitten, -1, died),
    individual::VariableUpdate$new(human, last_infected, -1, died),
    individual::VariableUpdate$new(human, infection_schedule, -1, died),
    individual::VariableUpdate$new(human, asymptomatic_infection_schedule, -1, died),
    individual::VariableUpdate$new(human, icm, birth_icm, died),
    individual::VariableUpdate$new(human, ivm, birth_ivm, died),
    individual::VariableUpdate$new(human, ib, -1, died),
    individual::VariableUpdate$new(human, ica, -1, died),
    individual::VariableUpdate$new(human, iva, -1, died),
    individual::VariableUpdate$new(human, id, -1, died),
    individual::VariableUpdate$new(human, is_severe, 0, died)
  )
}

died_from_severe <- function(severe, untreated, treated, ftv, v) {
  unsuccessful_treatment <- treated[bernoulli(length(treated), ftv)]
  at_risk <- intersect(
    severe,
    c(untreated, treated[unsuccessful_treatment])
  )
  at_risk[bernoulli(length(at_risk), v)]
}

died_naturally <- function(age, table) {
  age[age > 100] <- 100
  probability <- table[age]
  which(bernoulli(length(age), probability))
}

sample_mothers <- function(sampleable, died, groups) {
  vnapply(
    groups[died],
    function(source_group) {
      sample(which(sampleable & (groups == source_group)), 1)[[1]]
    }
  )
}

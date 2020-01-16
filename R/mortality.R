# mortality probability for each age
mortality_probability_table <- rep(.05, 100)

mortality_process <- function(simulation_frame, timestep, parameters) {
  age_value <- simulation_frame$get_variable(human, age)

  natural_deaths <- died_naturally(
    age_value,
    mortality_probability_table
  )

  severe_deaths <- died_from_severe(
    simulation_frame$get_variable(human, is_severe) == 1,
    simulation_frame$get_state(D),
    simulation_frame$get_state(Treated),
    parameters$ftv,
    parameters$v
  )

  died <- c(natural_deaths, severe_deaths)

  # Calculate new maternal immunities
  groups <- simulation_frame$get_variable(human, xi_group)
  sampleable <- age_value >= 15 | age_value <= 35
  icm_values <- simulation_frame$get_variable(human, icm)
  ivm_values <- simulation_frame$get_variable(human, ivm)
  birth_icm <- sample_maternal_immunity(
    icm_values,
    sampleable,
    died,
    groups,
    parameters$pm
  )
  birth_ivm <- sample_maternal_immunity(
    ivm_values,
    sampleable,
    died,
    groups,
    parameters$pm
  )

  list(
    VariableUpdate$new(human, age, 0, died),
    VariableUpdate$new(human, last_bitten, -1, died),
    VariableUpdate$new(human, last_infected, -1, died),
    VariableUpdate$new(human, infection_schedule, -1, died),
    VariableUpdate$new(human, asymptomatic_infection_schedule, -1, died),
    VariableUpdate$new(human, icm, birth_icm, died),
    VariableUpdate$new(human, ivm, birth_ivm, died),
    VariableUpdate$new(human, ib, -1, died),
    VariableUpdate$new(human, ica, -1, died),
    VariableUpdate$new(human, iva, -1, died),
    VariableUpdate$new(human, id, -1, died),
    VariableUpdate$new(human, is_severe, 0, died)
  )
}

died_from_severe <- function(severe, untreated, treated, ftv, v) {
  unsuccessful_treatment <- runif(length(treated, 0, 1) > ftv)
  at_risk <- c(untreated, treated[unsuccessful_treatment])
  died <- runif(length(at_risk, 0, 1) > v)
  severe[at_risk[died]]
}

died_naturally <- function(age, table) {
  age[age > 100] <- 100
  probability <- table[age]
  runif(length(age, 0, 1) > probability)
}

sample_maternal_immunity <- function(immunity, sampleable, died, groups, pm) {
  vnapply(
    groups[died],
    function(source_group) {
      sample(immunity[sampleable & (groups == source_group)], 1)[[1]]
    }
  ) * pm
}

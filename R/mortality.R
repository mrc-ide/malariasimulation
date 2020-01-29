#' @title Human mortality process
#' @description
#' This is the process for human mortality, it defines which humans die from
#' natural causes and severe infection and replaces dead individuals with
#' newborns.
create_mortality_process <- function(human, D, Treated, variables) {
  function(simulation_frame, timestep, parameters) {
    age <- simulation_frame$get_variable(human, variables$age)

    natural_deaths <- died_naturally(
      age,
      parameters$mortality_probability_table
    )

    severe_deaths <- died_from_severe(
      which(simulation_frame$get_variable(human, variables$is_severe) == 1),
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
    groups <- simulation_frame$get_variable(human, variables$xi_group)
    sampleable <- age >= 15 | age <= 35
    icm <- simulation_frame$get_variable(human, variables$icm)
    ivm <- simulation_frame$get_variable(human, variables$ivm)
    mothers <- sample_mothers(sampleable, died, groups)
    birth_icm <- icm[mothers] * parameters$pm
    birth_ivm <- ivm[mothers] * parameters$pm

    list(
      individual::VariableUpdate$new(human, variables$age, 0, died),
      individual::VariableUpdate$new(human, variables$last_bitten, -1, died),
      individual::VariableUpdate$new(human, variables$last_infected, -1, died),
      individual::VariableUpdate$new(human, variables$infection_schedule, -1, died),
      individual::VariableUpdate$new(human, variables$asymptomatic_infection_schedule, -1, died),
      individual::VariableUpdate$new(human, variables$icm, birth_icm, died),
      individual::VariableUpdate$new(human, variables$ivm, birth_ivm, died),
      individual::VariableUpdate$new(human, variables$ib, -1, died),
      individual::VariableUpdate$new(human, variables$ica, -1, died),
      individual::VariableUpdate$new(human, variables$iva, -1, died),
      individual::VariableUpdate$new(human, variables$id, -1, died),
      individual::VariableUpdate$new(human, variables$is_severe, 0, died)
    )
  }
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

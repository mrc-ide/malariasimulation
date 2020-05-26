#' @title Human mortality process
#' @description
#' This is the process for human mortality, it defines which humans die from
#' natural causes and severe infection and replaces dead individuals with
#' newborns.
#' @param human, the human individual
#' @param D, the diseased state
#' @param variables, the model variables to reset
create_mortality_process <- function(human, D, variables, events) {
  function(api) {
    parameters <- api$get_parameters()
    age <- api$get_variable(human, variables$age)

    # If this is the first timestep, we have to kick off the aging event
    if (api$get_timestep() == 1) {
      api$schedule(events$birthday, seq_len(parameters$human_population), 365)
    }

    natural_deaths <- died_naturally(
      age,
      parameters$mortality_probability_table
    )

    severe_deaths <- died_from_severe(
      which(api$get_variable(human, variables$is_severe) == 1),
      api$get_state(human, D),
      parameters$v
    )

    died <- c(natural_deaths, severe_deaths)

    # Calculate new maternal immunities
    groups <- api$get_variable(human, variables$xi_group)
    sampleable <- age >= 15 | age <= 35
    icm <- api$get_variable(human, variables$icm)
    ivm <- api$get_variable(human, variables$ivm)
    mothers <- sample_mothers(sampleable, died, groups)
    birth_icm <- icm[mothers] * parameters$pcm
    birth_ivm <- ivm[mothers] * parameters$pvm

    api$clear_schedule(events$infection, died)
    api$clear_schedule(events$asymptomatic_infection, died)
    api$clear_schedule(events$subpatent_infection, died)
    api$clear_schedule(events$recovery, died)
    api$clear_schedule(events$birthday, died)
    api$schedule(events$birthday, died, 365)

    api$queue_variable_update(human, variables$age, 0, died)
    api$queue_variable_update(human, variables$last_bitten, -1, died)
    api$queue_variable_update(human, variables$last_infected, -1, died)
    api$queue_variable_update(human, variables$icm, birth_icm, died)
    api$queue_variable_update(human, variables$ivm, birth_ivm, died)
    api$queue_variable_update(human, variables$ib, 0, died)
    api$queue_variable_update(human, variables$ica, 0, died)
    api$queue_variable_update(human, variables$iva, 0, died)
    api$queue_variable_update(human, variables$id, 0, died)
    api$queue_variable_update(human, variables$is_severe, 0, died)
  }
}

died_from_severe <- function(severe, diseased, v) {
  at_risk <- intersect(severe, diseased)
  at_risk[bernoulli(length(at_risk), v)]
}

died_naturally <- function(age, table) {
  age[age > 99] <- 99
  probability <- table[age + 1]
  which(bernoulli(length(age), probability))
}

sample_mothers <- function(sampleable, died, groups) {
  group_lookup <- list()
  for (group in unique(groups)) {
    group_lookup[[group]] <- which(sampleable & (groups == group))
  }
  vnapply(
    groups[died],
    function(source_group) {
      sample(group_lookup[[source_group]], 1)[[1]]
    }
  )
}

#' @title Human mortality process
#' @description
#' This is the process for human mortality, it defines which humans die from
#' natural causes and severe infection and replaces dead individuals with
#' newborns.
#' @param human the human individual
#' @param D the diseased state
#' @param variables the model variables to reset
#' @param events the model events to reset
create_mortality_process <- function(human, D, variables, events) {
  function(api) {
    parameters <- api$get_parameters()
    timestep <- api$get_timestep()
    age <- get_age(
      api$get_variable(human, variables$birth),
      timestep
    ) / 365

    natural_deaths <- which(
      bernoulli(
        parameters$human_population, parameters$mortality_rate
      )
    )
    
    api$render('natural_deaths', length(natural_deaths))

    severe_deaths <- died_from_severe(
      which(api$get_variable(human, variables$is_severe) == 1),
      api$get_state(human, D),
      parameters$v
    )

    api$render('severe_deaths', length(severe_deaths))

    died <- union(natural_deaths, severe_deaths)

    if (length(died) > 0) {
      # Calculate new maternal immunities
      groups <- api$get_variable(human, variables$xi_group)
      sampleable <- age >= 15 | age <= 35
      ica <- api$get_variable(human, variables$ica)
      iva <- api$get_variable(human, variables$iva)
      mothers <- sample_mothers(sampleable, died, groups)
      birth_icm <- ica[mothers] * parameters$pcm
      birth_ivm <- iva[mothers] * parameters$pvm

      api$clear_schedule(events$infection, died)
      api$clear_schedule(events$asymptomatic_infection, died)

      api$queue_variable_update(human, variables$birth, timestep, died)
      api$queue_variable_update(human, variables$last_bitten, -1, died)
      api$queue_variable_update(human, variables$last_infected, -1, died)
      api$queue_variable_update(human, variables$ib, 0, died)
      api$queue_variable_update(human, variables$ica, 0, died)
      api$queue_variable_update(human, variables$iva, 0, died)
      api$queue_variable_update(human, variables$id, 0, died)
      api$queue_variable_update(human, variables$icm, birth_icm, died)
      api$queue_variable_update(human, variables$ivm, birth_ivm, died)
      # xi and xi group survive rebirth
    }
  }
}

died_from_severe <- function(severe, diseased, v) {
  at_risk <- intersect(severe, diseased)
  at_risk[bernoulli(length(at_risk), v)]
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

#' @title Human mortality process
#' @description
#' This is the process for human mortality, it defines which humans die from
#' natural causes and severe infection and replaces dead individuals with
#' newborns.
#' @param human the human individual
#' @param D the diseased state
#' @param Tr the treated state
#' @param variables the model variables to reset
#' @param events the model events to reset
create_mortality_process <- function(human, D, Tr, variables, events) {
  function(api) {
    parameters <- api$get_parameters()
    timestep <- api$get_timestep()
    age <- get_age(
      api$get_variable(human, variables$birth),
      timestep
    ) / 365

    died <- which(bernoulli(
        parameters$human_population, parameters$mortality_rate
    ))
    
    api$render('natural_deaths', length(died))

    if (parameters$severe_enabled == 1) {
      severe_deaths <- died_from_severe(
        which(api$get_variable(human, variables$is_severe) == 1),
        api$get_state(human, D),
        parameters$v,
        api$get_state(human, Tr),
        parameters$fvt
      )
      api$render('severe_deaths', length(severe_deaths))
      died <- union(died, severe_deaths)
    }

    if (length(died) > 0) {
      # Calculate new maternal immunities
      groups <- api$get_variable(human, variables$zeta_group)
      sampleable <- age >= 15 | age <= 35
      ica <- api$get_variable(human, variables$ica)
      iva <- api$get_variable(human, variables$iva)
      mothers <- sample_mothers(sampleable, died, groups)
      birth_icm <- ica[mothers] * parameters$pcm
      birth_ivm <- iva[mothers] * parameters$pvm

      api$clear_schedule(events$infection, died)
      api$clear_schedule(events$asymptomatic_infection, died)

      api$queue_variable_update(human, variables$birth, timestep, died)
      api$queue_variable_update(human, variables$last_boosted_ib, -1, died)
      api$queue_variable_update(human, variables$last_boosted_ica, -1, died)
      api$queue_variable_update(human, variables$last_boosted_iva, -1, died)
      api$queue_variable_update(human, variables$last_boosted_id, -1, died)
      api$queue_variable_update(human, variables$ib, 0, died)
      api$queue_variable_update(human, variables$ica, 0, died)
      api$queue_variable_update(human, variables$iva, 0, died)
      api$queue_variable_update(human, variables$id, 0, died)
      api$queue_variable_update(human, variables$icm, birth_icm, died)
      api$queue_variable_update(human, variables$ivm, birth_ivm, died)
      api$queue_variable_update(human, variables$drug, 0, died)
      api$queue_variable_update(human, variables$drug_time, -1, died)
      api$queue_variable_update(human, variables$infectivity, 0, died)
      # zeta and zeta group survive rebirth
    }
  }
}

died_from_severe <- function(severe, diseased, v, treated, fvt) {
  at_risk <- intersect(severe, diseased)
  treated_at_risk <- intersect(severe, treated)
  c(
    at_risk[bernoulli(length(at_risk), v)],
    treated_at_risk[bernoulli(length(treated_at_risk), fvt)]
  )
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

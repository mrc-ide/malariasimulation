#' @title Human mortality process
#' @description
#' This is the process for human mortality, it defines which humans die from
#' natural causes and severe infection and replaces dead individuals with
#' newborns.
#' @param variables the model variables to reset
#' @param events the model events to reset
#' @param renderer the model renderer
#' @param parameters model parameters
create_mortality_process <- function(variables, events, renderer, parameters) {
  function(timestep) {
    age <- get_age(
      variables$birth$get_values(),
      timestep
    ) / 365

    died <- bernoulli(
      parameters$human_population, parameters$mortality_rate
    )
    
    renderer$render('natural_deaths', length(died))

    if (parameters$severe_enabled == 1) {
      severe_deaths <- died_from_severe(
        which(variables$is_severe$get_values() == 1),
        variables$state$get_index_of('D'),
        parameters$v,
        variables$state$get_index_of('Tr'),
        parameters$fvt
      )
      renderer$render('severe_deaths', length(severe_deaths))
      died <- union(died, severe_deaths)
    }

    if (length(died) > 0) {
      # Calculate new maternal immunities
      groups <- variables$zeta_group$get_values()
      sampleable <- age >= 15 | age <= 35
      ica <- variables$ica$get_values()
      iva <- variables$iva$get_values()
      mothers <- sample_mothers(sampleable, died, groups)
      birth_icm <- ica[mothers] * parameters$pcm
      birth_ivm <- iva[mothers] * parameters$pvm

      events$infection$clear_schedule(died)
      events$asymptomatic_infection$clear_schedule(died)

      variables$birth$queue_update(timestep, died)
      variables$last_boosted_ib$queue_update(-1, died)
      variables$last_boosted_ica$queue_update(-1, died)
      variables$last_boosted_iva$queue_update(-1, died)
      variables$last_boosted_id$queue_update(-1, died)
      variables$ib$queue_update(0, died)
      variables$ica$queue_update(0, died)
      variables$iva$queue_update(0, died)
      variables$id$queue_update(0, died)
      variables$icm$queue_update(birth_icm, died)
      variables$ivm$queue_update(birth_ivm, died)
      variables$drug$queue_update(0, died)
      variables$drug_time$queue_update(-1, died)
      variables$infectivity$queue_update(0, died)
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

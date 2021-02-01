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

    died <- individual::Bitset$new(parameters$human_population)
    died <- sample_bitset(died$not(), parameters$mortality_rate)

    renderer$render('natural_deaths', died$size(), timestep)

    if (parameters$severe_enabled == 1) {
      severe_deaths <- died_from_severe(
        variables$is_severe$get_index_of('yes'),
        variables$state$get_index_of('D'),
        parameters$v,
        variables$state$get_index_of('Tr'),
        parameters$fvt
      )
      renderer$render('severe_deaths', severe_deaths$size(), timestep)
      died$or(severe_deaths)
    }

    if (died$size() > 0) {
      # Calculate new maternal immunities
      sampleable <- individual::Bitset$new(parameters$human_population)
      sampleable$insert(which(age >= 15 | age <= 35))
      mothers <- sample_mothers(
        sampleable,
        died,
        parameters$n_heterogeneity_groups,
        variables$zeta_group
      )
      ica <- variables$ica$get_values()
      iva <- variables$iva$get_values()
      birth_icm <- ica[mothers] * parameters$pcm
      birth_ivm <- iva[mothers] * parameters$pvm

      events$infection$clear_schedule(died)
      events$asymptomatic_infection$clear_schedule(died)

      variables$state$queue_update('S', died)
      variables$is_severe$queue_update('no', died)
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
  at_risk <- severe$copy()$and(diseased)
  treated_at_risk <- severe$copy()$and(treated)
  sample_bitset(at_risk, v)$or(sample_bitset(treated_at_risk, fvt))
}

sample_mothers <- function(sampleable, died, n, groups) {
  group_lookup <- list()
  for (group in seq(n)) {
    group_lookup[[group]] <- sampleable$and(
      groups$get_index_of(toString(group))
    )$to_vector()
  }
  vnapply(
    groups[died],
    function(source_group) {
      sample(group_lookup[[source_group]], 1)[[1]]
    }
  )
}

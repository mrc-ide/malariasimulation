#' @title Human mortality process
#' @description
#' This is the process for human mortality, it defines which humans die from
#' natural causes and severe infection and replaces dead individuals with
#' newborns.
#' @param variables the model variables to reset
#' @param events the model events to reset
#' @param renderer the model renderer
#' @param parameters model parameters
#' @noRd
create_mortality_process <- function(variables, events, renderer, parameters) {
  function(timestep) {
    age <- trunc(get_age(
      variables$birth$get_values(),
      timestep
    ) / 365)

    died <- individual::Bitset$new(parameters$human_population)
    died <- sample_bitset(died$not(), 1 / parameters$average_age)

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
      # inherit immunity from parent in group
      sampleable <- individual::Bitset$new(parameters$human_population)
      sampleable$insert(which(age > 15 & age < 30))
      for (group in seq(parameters$n_heterogeneity_groups)) {

        # get the individuals who died in this group
        group_index <- variables$zeta_group$get_index_of(toString(group))
        died_group <- group_index$copy()$and(died)

        if (died_group$size() > 0) {
          # find their mothers
          potential_mothers <- group_index$and(sampleable)$to_vector()
          if (length(potential_mothers) == 0) {
            potential_mothers = seq(parameters$human_population)
          }
          mothers <- potential_mothers[
            sample.int(
              length(potential_mothers),
              died_group$size(),
              replace = TRUE
            )
          ]

          # set their maternal immunities
          birth_icm <- variables$ica$get_values(mothers) * parameters$pcm
          birth_ivm <- variables$ica$get_values(mothers) * parameters$pvm
          variables$icm$queue_update(birth_icm, died_group)
          variables$ivm$queue_update(birth_ivm, died_group)
        }
      }

      events$detection$clear_schedule(died)
      events$clinical_infection$clear_schedule(died)
      events$asymptomatic_infection$clear_schedule(died)
      events$asymptomatic_progression$clear_schedule(died)
      events$subpatent_progression$clear_schedule(died)
      events$recovery$clear_schedule(died)
      events$throw_away_net$clear_schedule(died)

      # new birthday
      variables$birth$queue_update(timestep, died)

      # non-maternal immunity
      variables$last_boosted_ib$queue_update(-1, died)
      variables$last_boosted_ica$queue_update(-1, died)
      variables$last_boosted_iva$queue_update(-1, died)
      variables$last_boosted_id$queue_update(-1, died)
      variables$ib$queue_update(0, died)
      variables$ica$queue_update(0, died)
      variables$iva$queue_update(0, died)
      variables$id$queue_update(0, died)
      variables$state$queue_update('S', died)

      # treatment
      variables$drug$queue_update(0, died)
      variables$drug_time$queue_update(-1, died)

      # vaccination
      variables$rtss_vaccinated$queue_update(-1, died)
      variables$rtss_boosted$queue_update(-1, died)
      variables$tbv_vaccinated$queue_update(-1, died)

      # onwards infectiousness
      variables$infectivity$queue_update(0, died)

      # vector control
      variables$net_time$queue_update(-1, died)
      variables$spray_time$queue_update(-1, died)

      # misc.
      variables$is_severe$queue_update('no', died)
      # zeta and zeta group survive rebirth
    }
  }
}

died_from_severe <- function(severe, diseased, v, treated, fvt) {
  at_risk <- severe$copy()$and(diseased)
  treated_at_risk <- severe$copy()$and(treated)
  sample_bitset(at_risk, v)$or(sample_bitset(treated_at_risk, fvt))
}

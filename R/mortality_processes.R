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
    at_risk <- variables$state$get_index_of('NonExistent')$not(TRUE)

    if (is.null(parameters$deathrate_timesteps)) {
      died <- sample_bitset(at_risk, 1 / parameters$average_age)
      renderer$render('natural_deaths', died$size(), timestep)
    } else {
      age <- get_age(variables$birth$get_values(at_risk), timestep)
      last_deathrate <- which(timestep >= parameters$deathrate_timesteps)[[-1]]
      age_groups <- .bincode(age, c(0, parameters$deathrate_agegroups))
      deathrates <- parameters$deathrates[age_groups, last_deathrate]
      deathrates[is.na(deathrates)] <- 1
      died <- bitset_at(at_risk, bernoulli_multi_p(deathrates))
      renderer$render('natural_deaths', died$size(), timestep)
    }
    reset_target(variables, events, died, 'NonExistent')
  }
}

create_birth_process <- function(variables, events, parameters) {
  function(timestep) {
    at_risk <- variables$state$get_index_of('NonExistent')
    if (is.null(parameters$birthrate_timesteps)) {
      born <- sample_bitset(at_risk, 1 / parameters$average_age)
    } else {
      last_birthrate <- which(timestep >= parameters$birthrate_timesteps)[[-1]]
      born <- sample_bitset(at_risk, parameters$birthrates[last_birthrate])
    }
    sample_maternal_immunity(variables, born, timestep)
    reset_target(variables, events, born, 'S')
  }
}

sample_maternal_immunity <- function(variables, target, timestep) {
  if (target$size() > 0) {
    age <- get_age(variables$birth$get_values(), timestep)
    # inherit immunity from parent in group
    sampleable <- individual::Bitset$new(parameters$human_population)
    sampleable$insert(which(trunc(age / 365) == 20))
    for (group in seq(parameters$n_heterogeneity_groups)) {

      # get the individuals who died in this group
      group_index <- variables$zeta_group$get_index_of(toString(group))
      target_group <- group_index$copy()$and(target)

      if (target_group$size() > 0) {
        # find their mothers
        potential_mothers <- group_index$and(sampleable)$to_vector()
        if (length(potential_mothers) == 0) {
          potential_mothers = seq(parameters$human_population)
        }
        mothers <- potential_mothers[
          sample.int(
            length(potential_mothers),
            target_group$size(),
            replace = TRUE
          )
        ]

        # set their maternal immunities
        birth_icm <- variables$ica$get_values(mothers) * parameters$pcm
        birth_ivm <- variables$ica$get_values(mothers) * parameters$pvm
        variables$icm$queue_update(birth_icm, target_group)
        variables$ivm$queue_update(birth_ivm, target_group)
      }
    }
  }
}

reset_target <- function(variables, events, target, state) {
  if (target$size() > 0) {
    # clear events
    to_clear <- c(
      'asymptomatic_progression',
      'subpatent_progression',
      'recovery',
      'clinical_infection',
      'asymptomatic_infection',
      'detection',
      'throw_away_net',
      'rtss_mass_doses',
      'rtss_mass_booster',
      'rtss_epi_doses',
      'rtss_epi_boosters'
    )
    for (event in unlist(events[to_clear])) {
      event$clear_schedule(target)
    }

    # new birthday
    variables$birth$queue_update(timestep, target)

    # non-maternal immunity
    variables$last_boosted_ib$queue_update(-1, target)
    variables$last_boosted_ica$queue_update(-1, target)
    variables$last_boosted_iva$queue_update(-1, target)
    variables$last_boosted_id$queue_update(-1, target)
    variables$ib$queue_update(0, target)
    variables$ica$queue_update(0, target)
    variables$iva$queue_update(0, target)
    variables$id$queue_update(0, target)
    variables$state$queue_update(state, target)

    # treatment
    variables$drug$queue_update(0, target)
    variables$drug_time$queue_update(-1, target)

    # vaccination
    variables$rtss_vaccinated$queue_update(-1, target)
    variables$rtss_boosted$queue_update(-1, target)
    variables$tbv_vaccinated$queue_update(-1, target)

    # onwards infectiousness
    variables$infectivity$queue_update(0, target)

    # vector control
    variables$net_time$queue_update(-1, target)
    variables$spray_time$queue_update(-1, target)
    # zeta and zeta group survive rebirth
  }
}

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

    pop <- get_human_population(parameters, timestep)
    if (!parameters$custom_demography) {
      died <- individual::Bitset$new(pop)$insert(
        bernoulli(pop, 1 / parameters$average_age)
      )
      renderer$render('natural_deaths', died$size(), timestep)
    } else {
      age <- get_age(variables$birth$get_values(), timestep)
      last_deathrate <- match_timestep(parameters$deathrate_timesteps, timestep)
      deathrates <- rep(1, pop)
      age_groups <- .bincode(age, c(0, parameters$deathrate_agegroups))
      in_range <- !is.na(age_groups)
      deathrates[in_range] <- parameters$deathrates[last_deathrate, age_groups[in_range]]
      died <- individual::Bitset$new(pop)$insert(bernoulli_multi_p(deathrates))
      renderer$render('natural_deaths', died$size(), timestep)
    }
    reset_target(variables, events, died, 'S', timestep)
    sample_maternal_immunity(variables, died, timestep, parameters)
  }
}

sample_maternal_immunity <- function(variables, target, timestep, parameters) {
  if (target$size() > 0) {
    pop <- get_human_population(parameters, timestep)
    age <- get_age(variables$birth$get_values(), timestep)
    # inherit immunity from parent in group
    sampleable <- individual::Bitset$new(pop)
    sampleable$insert(which(trunc(age / 365) == 20))
    for (group in seq(parameters$n_heterogeneity_groups)) {

      # get the individuals who died in this group
      group_index <- variables$zeta_group$get_index_of(toString(group))
      target_group <- group_index$copy()$and(target)

      if (target_group$size() > 0) {
        # find their mothers
        potential_mothers <- group_index$and(sampleable)$to_vector()
        if (length(potential_mothers) == 0) {
          mothers <- sample.int(
            pop,
            target_group$size(),
            replace = TRUE
          )
        } else {
          mothers <- potential_mothers[
            sample.int(
              length(potential_mothers),
              target_group$size(),
              replace = TRUE
            )
          ]
        }

        # set their maternal immunities
        birth_icm <- variables$ica$get_values(mothers) * parameters$pcm
        birth_ivm <- variables$ica$get_values(mothers) * parameters$pvm
        variables$icm$queue_update(birth_icm, target_group)
        variables$ivm$queue_update(birth_ivm, target_group)
      }
    }
  }
}

reset_target <- function(variables, events, target, state, timestep) {
  if (target$size() > 0) {
    # clear events
    to_clear <- c(
      'mass_pev_doses',
      'mass_pev_boosters',
      'pev_epi_doses',
      'pev_epi_boosters'
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
    variables$last_pev_timestep$queue_update(-1, target)
    variables$last_eff_pev_timestep$queue_update(-1, target)
    variables$pev_profile$queue_update(-1, target)
    variables$tbv_vaccinated$queue_update(-1, target)

    # onwards infectiousness
    variables$infectivity$queue_update(0, target)

    # zeta and zeta group and vector controls survive rebirth
  }
}

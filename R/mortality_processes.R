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
    reset_target(variables, events, died, 'S', parameters, timestep)
    sample_maternal_immunity(variables, died, timestep, parameters)
  }
}

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
create_verbose_mortality_process <- function(variables, events, renderer, parameters) {
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
    reset_target_verbose(variables, events, died, 'S', parameters, timestep)
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
        variables$icm$queue_update(birth_icm, target_group)
        if(parameters$parasite == "falciparum"){
          birth_ivm <- variables$iva$get_values(mothers) * parameters$pvm
          variables$ivm$queue_update(birth_ivm, target_group)
          
        } else if(parameters$parasite == "vivax"){
          birth_iam <- variables$iaa$get_values(mothers) * parameters$pcm
          variables$iam$queue_update(birth_iam, target_group)
        }
      }
    }
  }
}

reset_target <- function(variables, events, target, state, parameters, timestep) {
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
    variables$last_boosted_ica$queue_update(-1, target)
    variables$ica$queue_update(0, target)
    variables$state$queue_update(state, target)
    
    if(parameters$parasite == "falciparum"){
      variables$last_boosted_ib$queue_update(-1, target)
      variables$last_boosted_iva$queue_update(-1, target)
      variables$last_boosted_id$queue_update(-1, target)
      variables$ib$queue_update(0, target)
      variables$iva$queue_update(0, target)
      variables$id$queue_update(0, target)
      
    } else if (parameters$parasite == "vivax"){
      variables$last_boosted_iaa$queue_update(-1, target)
      variables$iaa$queue_update(0, target)
      variables$hypnozoites$queue_update(0, target)
    }

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
    variables$progression_rates$queue_update(0, target)
    
    # zeta and zeta group and vector controls survive rebirth
  }
}
reset_target_verbose <- function(variables, events, target, state, parameters, timestep) {
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
    
    if(parameters$mortality_verbose){
      # recording_people <- target$or(variables$birth$get_index_of(a = parameters$lower_age_bound, b = parameters$upper_age_bound))
      # states <- variables$state$get_values(recording_people$to_vector())
      # personal_inds <- variables$personal_tracker_index$get_values(recording_people$to_vector())
      states <- variables$state$get_values(target$to_vector())
      personal_inds <- variables$personal_tracker_index$get_values(target$to_vector())
      print_to_csv(parameters$file_name, timestep, personal_inds, "died", states, parameters$start_time)
    }

    quantity_to_update <- target$size()
    curr_max_ind <- max(variables$personal_tracker_index$get_values())
    variables$personal_tracker_index$queue_update(seq(curr_max_ind + 1, curr_max_ind + quantity_to_update), target)

    # new birthday
    variables$birth$queue_update(timestep, target)

    # non-maternal immunity
    variables$last_boosted_ica$queue_update(-1, target)
    variables$ica$queue_update(0, target)
    variables$state$queue_update(state, target)
    
    if(parameters$parasite == "falciparum"){
      variables$last_boosted_ib$queue_update(-1, target)
      variables$last_boosted_iva$queue_update(-1, target)
      variables$last_boosted_id$queue_update(-1, target)
      variables$ib$queue_update(0, target)
      variables$iva$queue_update(0, target)
      variables$id$queue_update(0, target)
      
    } else if (parameters$parasite == "vivax"){
      variables$last_boosted_iaa$queue_update(-1, target)
      variables$iaa$queue_update(0, target)
      variables$hypnozoites$queue_update(0, target)
    }

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
    variables$progression_rates$queue_update(0, target)
    
    if(parameters$mortality_verbose){
      states <- variables$state$get_values(target$to_vector())
      print_to_csv(parameters$file_name, timestep + 1, seq(curr_max_ind + 1, curr_max_ind + quantity_to_update), "born", states, parameters$start_time)
    }
    # zeta and zeta group and vector controls survive rebirth
  }
}

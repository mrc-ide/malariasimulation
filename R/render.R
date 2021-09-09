in_age_range <- function(birth, timestep, lower, upper) {
  birth$get_index_of(a = timestep - upper, b = timestep - lower)
}

#' @title Render prevalence statistics
#' 
#' @description renders prevalence numerators and denominators for indivduals
#' detected by microscopy and with severe malaria
#' 
#' @param state human infection state
#' @param birth variable for birth of the individual
#' @param is_severe variable for if individual has severe malaria
#' @param immunity to detection
#' @param parameters model parameters
#' @param renderer model renderer
#' 
#' @noRd
create_prevelance_renderer <- function(
  state,
  birth,
  is_severe,
  immunity,
  parameters,
  renderer
  ) {
  function(timestep) {
    age <- get_age(birth$get_values(), timestep)
    asymptomatic <- state$get_index_of('A')
    asymptomatic_vector <- asymptomatic$to_vector()
    asymptomatic_detected <- bernoulli_multi_p(
      probability_of_detection(
        age[asymptomatic_vector],
        immunity$get_values(asymptomatic),
        parameters
      )
    )
    detected <- state$get_index_of(c('Tr', 'D'))$or(
      bitset_at(asymptomatic, asymptomatic_detected)
    )

    severe <- is_severe$get_index_of('yes')$and(detected)
    for (i in seq_along(parameters$prevalence_rendering_min_ages)) {
      lower <- parameters$prevalence_rendering_min_ages[[i]]
      upper <- parameters$prevalence_rendering_max_ages[[i]]
      in_age <- in_age_range(birth, timestep, lower, upper)
      renderer$render(
        paste0('n_', lower, '_', upper),
        in_age$size(),
        timestep
      ) 
      renderer$render(
        paste0('n_detect_', lower, '_', upper),
        in_age$and(detected)$size(),
        timestep
      )
    }
    for (i in seq_along(parameters$severe_prevalence_rendering_min_ages)) {
      lower <- parameters$severe_prevalence_rendering_min_ages[[i]]
      upper <- parameters$severe_prevalence_rendering_max_ages[[i]]
      in_age <- in_age_range(birth, timestep, lower, upper)
      renderer$render(
        paste0('n_', lower, '_', upper),
        in_age$size(),
        timestep
      ) 
      renderer$render(
        paste0('n_severe_', lower, '_', upper),
        in_age$and(severe)$size(),
        timestep
      )
    }
  }
}

#' @title Render incidence statistics
#' 
#' @description renders incidence (new for this timestep) for indivduals
#' detected by microscopy and with severe malaria
#' 
#' @param birth variable for birth of the individual
#' @param is_severe variable for if individual has severe malaria
#' @param parameters model parameters
#' @param renderer model renderer
#' @param timestep current target
#' @param target newly infected population
#' 
#' @noRd
incidence_renderer <- function(
  birth,
  is_severe,
  parameters,
  renderer,
  timestep,
  target
  ) {
  severe <- is_severe$get_index_of('yes')$and(target)
  for (i in seq_along(parameters$incidence_rendering_min_ages)) {
    lower <- parameters$incidence_rendering_min_ages[[i]]
    upper <- parameters$incidence_rendering_max_ages[[i]]
    in_age <- in_age_range(birth, timestep, lower, upper)
    renderer$render(paste0('n_', lower, '_', upper), in_age$size(), timestep)
    renderer$render(
      paste0('n_inc_', lower, '_', upper),
      in_age$and(target)$size(),
      timestep
    )
  }
  for (i in seq_along(parameters$severe_incidence_rendering_min_ages)) {
    lower <- parameters$severe_incidence_rendering_min_ages[[i]]
    upper <- parameters$severe_incidence_rendering_max_ages[[i]]
    in_age <- in_age_range(birth, timestep, lower, upper)
    renderer$render(paste0('n_', lower, '_', upper), in_age$size(), timestep)
    renderer$render(
      paste0('n_inc_severe_', lower, '_', upper),
      in_age$and(target)$and(severe)$size(),
      timestep
    )
  }
}

#' @title Render clinical incidence statistics
#' 
#' @description renders clinical incidence (new for this timestep)
#' 
#' @param birth variable for birth of the individual
#' @param parameters model parameters
#' @param renderer model renderer
#' 
#' @noRd
clinical_incidence_renderer <- function(
  birth,
  parameters,
  renderer,
  target,
  timestep
  ) {
  for (i in seq_along(parameters$clinical_incidence_rendering_min_ages)) {
    lower <- parameters$clinical_incidence_rendering_min_ages[[i]]
    upper <- parameters$clinical_incidence_rendering_max_ages[[i]]
    in_age <- in_age_range(birth, timestep, lower, upper)
    renderer$render(paste0('n_', lower, '_', upper), in_age$size(), timestep)
    renderer$render(
      paste0('n_inc_clinical_', lower, '_', upper),
      in_age$and(target)$size(),
      timestep
    )
  }
}

create_variable_mean_renderer_process <- function(
  renderer,
  names,
  variables
) {
  function(timestep) {
    for (i in seq_along(variables)) {
      renderer$render(
        paste0(names[[i]], '_mean'),
        mean(variables[[i]]$get_values()),
        timestep
      )
    }
  }
}

create_total_M_renderer_individual <- function(
  mosquito_state,
  species,
  renderer,
  parameters
  ) {
  function(timestep) {
    adult <- mosquito_state$get_index_of('NonExistent')$not()
    for (i in seq_along(parameters$species)) {
      renderer$render(
        paste0('total_M_', i),
        species$get_index_of(parameters$species[[i]])$and(adult)$size(),
        timestep
      )
    }

    renderer$render('total_M', adult$size(), timestep)
  }
}

create_total_M_renderer_compartmental <- function(renderer, solvers) {
  function(timestep) {
    total_M <- 0
    for (i in seq_along(solvers)) {
      row <- solver_get_states(solvers[[i]])
      species_M <- sum(row[ADULT_ODE_INDICES])
      total_M <- total_M + species_M
      renderer$render(paste0('total_M_', i), species_M, timestep)
    }

    renderer$render('total_M', total_M, timestep)
  }
}

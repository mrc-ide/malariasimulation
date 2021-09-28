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
#' @param renderer object for model outputs
#' @param target incidence population
#' @param prefix for model outputs
#' @param lowers age bounds
#' @param uppers age bounds
#' @param timestep current target
#' 
#' @noRd
incidence_renderer <- function(
  birth,
  renderer,
  target,
  prefix,
  lowers,
  uppers,
  timestep
  ) {
  for (i in seq_along(lowers)) {
    lower <- lowers[[i]]
    upper <- uppers[[i]]
    in_age <- in_age_range(birth, timestep, lower, upper)
    renderer$render(paste0('n_', lower, '_', upper), in_age$size(), timestep)
    renderer$render(
      paste0(prefix, lower, '_', upper),
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

create_vector_count_renderer_individual <- function(
  mosquito_state,
  species,
  state,
  renderer,
  parameters
  ) {
  function(timestep) {
    adult <- mosquito_state$get_index_of('NonExistent')$not(TRUE)
    for (i in seq_along(parameters$species)) {
      species_name <- parameters$species[[i]]
      species_index <- species$get_index_of(species_name)
      for (s in state$get_categories()) {
        renderer$render(
          paste0(s, '_', species_name, '_count'),
          state$get_index_of(s)$and(species_index)$size(),
          timestep
        )
      }
    }
  }
}

create_total_M_renderer_compartmental <- function(renderer, solvers, parameters) {
  function(timestep) {
    total_M <- 0
    for (i in seq_along(solvers)) {
      row <- solver_get_states(solvers[[i]])
      species_M <- sum(row[ADULT_ODE_INDICES])
      total_M <- total_M + species_M
      renderer$render(paste0('total_M_', parameters$species[[i]]), species_M, timestep)
    }
  }
}

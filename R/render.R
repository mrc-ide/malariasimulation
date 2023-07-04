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
#' @param immunity to detection
#' @param parameters model parameters
#' @param renderer model renderer
#' 
#' @noRd
create_prevelance_renderer <- function(
  state,
  birth,
  immunity,
  parameters,
  renderer
  ) {
  function(timestep) {
    
    asymptomatic <- state$get_index_of('A')
    
    if(parameters$parasite =="falciparum"){
      prob <- probability_of_detection(
        get_age(birth$get_values(asymptomatic), timestep),
        immunity$get_values(asymptomatic),
        parameters
      )
      asymptomatic_detected <- bitset_at(asymptomatic, bernoulli_multi_p(prob))
      
    } else if (parameters$parasite =="vivax"){
      # The vivax model defines asymptomatic infections as being detectable by
      # light microscopy
      prob <- rep(1, asymptomatic$size())
      asymptomatic_detected <- state$get_index_of('A')
    }
    
    clinically_detected <- state$get_index_of(c('Tr', 'D'))
    detected <- clinically_detected$copy()$or(asymptomatic_detected)

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
        in_age$copy()$and(detected)$size(),
        timestep
      )
      renderer$render(
        paste0('p_detect_', lower, '_', upper),
        in_age$copy()$and(clinically_detected)$size() + sum(
          prob[bitset_index(asymptomatic, in_age)]
        ),
        timestep
      )
    }
  }
}

#' @title Render incidence statistics
#' 
#' @description renders incidence (new for this timestep) for indivduals
#' 
#' @param birth variable for birth of the individual
#' @param renderer object for model outputs
#' @param target incidence population
#' @param source_pop the population which is sampled for infection
#' @param prob probability of infection
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
  source_pop,
  prob,
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
      paste0('n_', prefix, lower, '_', upper),
      in_age$copy()$and(target)$size(),
      timestep
    )

    renderer$render(
      paste0('p_', prefix, lower, '_', upper),
      sum(prob[bitset_index(source_pop, in_age)]),
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

create_age_group_renderer <- function(
    birth,
    parameters,
    renderer
) {
  function(timestep) {
    for (i in seq_along(parameters$age_group_rendering_min_ages)) {
      lower <- parameters$age_group_rendering_min_ages[[i]]
      upper <- parameters$age_group_rendering_max_ages[[i]]
      in_age <- in_age_range(birth, timestep, lower, upper)
      renderer$render(
        paste0('n_age_', lower, '_', upper),
        in_age$size(),
        timestep
      ) 
    }
  }
}

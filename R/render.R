count_of_subpopulation <- function(target, age, lower, upper) {
  in_range <- individual::Bitset$new(target$max_size)$insert(
    which((age >= lower) & (age <= upper))
  )
  n_in_range <- in_range$size()
  if (n_in_range == 0) {
    return(0)
  }
  n_in_range
}

count_of_detected <- function(target, age, lower, upper) {
  in_range <- individual::Bitset$new(target$max_size)$insert(
    which((age >= lower) & (age <= upper))
  )
  n_in_range <- in_range$size()
  if (n_in_range == 0) {
    return(0)
  }
  in_range$and(target)$size() 
}

# calculate # of people in age group, 
# number detected by microscopy (all), 
# number of severe (all)
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
      p <- count_of_subpopulation( 
        detected,
        age,
        lower,
        upper
      )
      renderer$render(paste0('n_', lower, '_', upper), p, timestep) 
      p <- count_of_detected( 
        detected,
        age,
        lower,
        upper
      )
      renderer$render(paste0('n_', 'detect_', lower, '_', upper), p, timestep) # NEW
    }
    for (i in seq_along(parameters$severe_prevalence_rendering_min_ages)) {
      lower <- parameters$severe_prevalence_rendering_min_ages[[i]]
      upper <- parameters$severe_prevalence_rendering_max_ages[[i]]
      p <- count_of_detected(
        severe,
        age,
        lower,
        upper
      )
      renderer$render(paste0('n_', 'severe_', lower, '_', upper), p, timestep)
    }
  }
}

# calculate number detected by microscopy (new in that timestep), 
# number of severe (new in that timestep)
create_incidence_renderer <- function(birth, is_severe, parameters, renderer) {
  function(timestep, target) {
    age <- get_age(birth$get_values(), timestep)
    severe <- is_severe$get_index_of('yes')$and(target)
    for (i in seq_along(parameters$incidence_rendering_min_ages)) {
      lower <- parameters$incidence_rendering_min_ages[[i]]
      upper <- parameters$incidence_rendering_max_ages[[i]]
      p <- count_of_detected( 
        target,
        age,
        lower,
        upper
      )
      renderer$render(paste0('n_', 'inc_', lower, '_', upper), p, timestep) 
    }
    for (i in seq_along(parameters$severe_incidence_rendering_min_ages)) {
      lower <- parameters$severe_incidence_rendering_min_ages[[i]]
      upper <- parameters$severe_incidence_rendering_max_ages[[i]]
      p <- count_of_detected(
        severe$copy()$and(target),
        age,
        lower,
        upper
      )
      renderer$render(paste0('n_', 'inc_', 'severe_', lower, '_', upper), p, timestep)
    }
  }
}

create_clinical_incidence_renderer <- function(birth, parameters, renderer) {
  function(timestep, target) {
    age <- get_age(birth$get_values(), timestep)
    for (i in seq_along(parameters$clinical_incidence_rendering_min_ages)) {
      lower <- parameters$clinical_incidence_rendering_min_ages[[i]]
      upper <- parameters$clinical_incidence_rendering_max_ages[[i]]
      p <- count_of_detected( 
        target,
        age,
        lower,
        upper
      )
      renderer$render(paste0('n_', 'inc_', 'clinical_', lower, '_', upper), p, timestep) 
    }
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

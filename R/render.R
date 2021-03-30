epi_from_subpopulation <- function(target, age, lower, upper) {
  in_range <- individual::Bitset$new(target$max_size)$insert(
    which((age >= lower) & (age <= upper))
  )

  n_in_range <- in_range$size()
  if (n_in_range == 0) {
    return(0)
  }
  in_range$and(target)$size() / n_in_range
}

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
      p <- epi_from_subpopulation(
        detected,
        age,
        lower,
        upper
      )
      renderer$render(paste0('pv_', lower, '_', upper), p, timestep)
    }
    for (i in seq_along(parameters$severe_prevalence_rendering_min_ages)) {
      lower <- parameters$severe_prevalence_rendering_min_ages[[i]]
      upper <- parameters$severe_prevalence_rendering_max_ages[[i]]
      p <- epi_from_subpopulation(
        severe,
        age,
        lower,
        upper
      )
      renderer$render(paste0('pv_severe_', lower, '_', upper), p, timestep)
    }
  }
}

create_incidence_renderer <- function(birth, is_severe, parameters, renderer) {
  function(timestep, target) {
    age <- get_age(birth$get_values(), timestep)
    severe <- is_severe$get_index_of('yes')$and(target)
    for (i in seq_along(parameters$incidence_rendering_min_ages)) {
      lower <- parameters$incidence_rendering_min_ages[[i]]
      upper <- parameters$incidence_rendering_max_ages[[i]]
      p <- epi_from_subpopulation(
        target,
        age,
        lower,
        upper
      )
      renderer$render(paste0('inc_', lower, '_', upper), p, timestep)
    }
    for (i in seq_along(parameters$severe_incidence_rendering_min_ages)) {
      lower <- parameters$severe_incidence_rendering_min_ages[[i]]
      upper <- parameters$severe_incidence_rendering_max_ages[[i]]
      p <- epi_from_subpopulation(
        severe,
        age,
        lower,
        upper
      )
      renderer$render(paste0('inc_severe_', lower, '_', upper), p, timestep)
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

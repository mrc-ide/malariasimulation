epi_from_subpopulation <- function(target, age, lower, upper) {
  in_range <- which((age >= lower) & (age <= upper))
  if (length(in_range) == 0) {
    return(0)
  }
  n_target <- length(intersect(target, in_range))
  n_target / length(in_range)
}

create_prevelance_renderer <- function(human, D, A, birth, is_severe) {
  function(api) {
    parameters <- api$get_parameters()
    infected <- api$get_state(human, D, A)
    severe <- intersect(
      infected,
      which(api$get_variable(human, is_severe) == 1)
    )
    age <- get_age(api$get_variable(human, birth), api$get_timestep())
    for (i in seq_along(parameters$prevalence_rendering_min_ages)) {
      lower <- parameters$prevalence_rendering_min_ages[[i]]
      upper <- parameters$prevalence_rendering_max_ages[[i]]
      p <- epi_from_subpopulation(
        infected,
        age,
        lower,
        upper
      )
      api$render(paste0('pv_', lower, '_', upper), p)
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
      api$render(paste0('pv_severe_', lower, '_', upper), p)
    }
  }
}

create_incidence_renderer <- function(human, birth, is_severe) {
  function(api, target) {
    parameters <- api$get_parameters()
    age <- get_age(api$get_variable(human, birth), api$get_timestep())
    severe <- intersect(
      target,
      which(api$get_variable(human, is_severe) == 1)
    )
    for (i in seq_along(parameters$incidence_rendering_min_ages)) {
      lower <- parameters$incidence_rendering_min_ages[[i]]
      upper <- parameters$incidence_rendering_max_ages[[i]]
      p <- epi_from_subpopulation(
        target,
        age,
        lower,
        upper
      )
      api$render(paste0('inc_', lower, '_', upper), p)
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
      api$render(paste0('inc_severe_', lower, '_', upper), p)
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

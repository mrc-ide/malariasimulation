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
    prob <- probability_of_detection(
      get_age(birth$get_values(asymptomatic), timestep),
      immunity$get_values(asymptomatic),
      parameters
    )
    asymptomatic_detected <- bitset_at(asymptomatic, bernoulli_multi_p(prob))

    clinically_detected <- state$get_index_of(c('Tr', 'D'))
    detected <- clinically_detected$copy()$or(asymptomatic_detected)

    for (i in seq_along(parameters$prevalence_rendering_min_ages)) {
      lower <- parameters$prevalence_rendering_min_ages[[i]]
      upper <- parameters$prevalence_rendering_max_ages[[i]]
      in_age <- in_age_range(birth, timestep, lower, upper)
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
#' @description renders incidence (new for this timestep) for individuals
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
    renderer$render(
      paste0('n_', prefix, lower, '_', upper),
      in_age$copy()$and(target)$size(),
      timestep
    )
  }
}

#' @title Render probability of incidence statistics
#' 
#' @description renders probability of incidence (new for this timestep) for individuals
#' 
#' @param birth variable for birth of the individual
#' @param renderer object for model outputs
#' @param source_pop the population which is sampled for infection (bitten SAU for incidence, infecte for clinical/severe)
#' @param prob probability of infection
#' @param prefix for model outputs
#' @param lowers age bounds
#' @param uppers age bounds
#' @param timestep current target
#' 
#' @noRd
incidence_probability_renderer <- function(
    birth,
    renderer,
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
      row <- solvers[[i]]$get_states()
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
    
    unique_age_combinations <- as.data.frame(unique(rbind(cbind(parameters$prevalence_rendering_min_ages, parameters$prevalence_rendering_max_ages),
          cbind(parameters$incidence_rendering_min_ages, parameters$incidence_rendering_max_ages),
          cbind(parameters$clinical_incidence_rendering_min_ages, parameters$clinical_incidence_rendering_max_ages),
          cbind(parameters$severe_incidence_rendering_min_ages, parameters$severe_incidence_rendering_max_ages),
          cbind(parameters$age_group_rendering_min_ages, parameters$age_group_rendering_max_ages))))
    
    ordered_unique_age_combinations <- unique_age_combinations[
      with(unique_age_combinations, order(V1, V2)),
    ]
    
    for (i in seq_along(ordered_unique_age_combinations$V1)) {
      lower <- ordered_unique_age_combinations$V1[[i]]
      upper <- ordered_unique_age_combinations$V2[[i]]
      in_age <- in_age_range(birth, timestep, lower, upper)
      renderer$render(
        paste0('n_', lower, '_', upper),
        in_age$size(),
        timestep
      ) 
    }
  }
}


populate_incidence_rendering_columns <- function(renderer, parameters){
  
  # infections must render in all simulations 
  renderer$set_default('n_infections', 0)
  
  # treatment associated only renders when drugs are used
  if(sum(unlist(parameters$clinical_treatment_coverages))>0){
    renderer$set_default('ft', 0)
    renderer$set_default('n_treated', 0)
    renderer$set_default('n_drug_efficacy_failures', 0)
    renderer$set_default('n_successfully_treated', 0)
  }
  
  # ETC, SPC only render when antimalarial resistance is on
  if(parameters$antimalarial_resistance){
    renderer$set_default('n_early_treatment_failure', 0)
    renderer$set_default('n_slow_parasite_clearance', 0)
  }
  
  if(length(parameters$incidence_rendering_min_ages)>0){
    render_initial_incidence(renderer,
                             parameters$incidence_rendering_min_ages,
                             parameters$incidence_rendering_max_ages,
                             "inc")
  }
  
  if(length(parameters$clinical_incidence_rendering_min_ages)>0){
    render_initial_incidence(renderer,
                             parameters$clinical_incidence_rendering_min_ages,
                             parameters$clinical_incidence_rendering_max_ages,
                             "inc_clinical")
  }

  if(length(parameters$severe_incidence_rendering_min_ages)>0){
    render_initial_incidence(renderer,
                             parameters$severe_incidence_rendering_min_ages,
                             parameters$severe_incidence_rendering_max_ages,
                             "inc_severe")
  }
  
}

render_initial_incidence <- function(renderer, lower_vals, upper_vals, inc_name){
  for (i in seq_along(lower_vals)){
    renderer$set_default(paste0('n_', inc_name, "_", lower_vals[i], '_', upper_vals[i]), 0)
    renderer$set_default(paste0('p_', inc_name, "_", lower_vals[i], '_', upper_vals[i]), 0)
  }
}
  
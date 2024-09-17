#' @title pmc process
#'
#' @description schedules individuals to be given perennial malaria 
#' chemoprevention according to an age-based strategy
#'
#' @param variables list of variables in the model
#' @param events a list of events in the model
#' @param parameters the model parameters
#' @param renderer a renderer object
#' @param correlations correlation parameters
#' @param coverages PMC coverage
#' @param timesteps PMC coverage change timesteps
#' @param drug PMC drug index
#' @noRd
create_pmc_process <- function(
    variables,
    events,
    parameters,
    renderer,
    correlations,
    coverages,
    timesteps,
    drug
){
  renderer$set_default('n_pmc_treated', 0)
  renderer$set_default(paste0('n_pmc_drug_efficacy_failures'), 0)
  renderer$set_default(paste0('n_pmc_successfully_treated'), 0)
  
  if(parameters$antimalarial_resistance){
    renderer$set_default(paste0('n_pmc_early_treatment_failure'), 0)
    renderer$set_default(paste0('n_pmc_slow_parasite_clearance'), 0)
  }
  
  function(timestep) {
    time_index <- match_timestep(ts = timesteps, t = timestep)
    if(time_index == 0){
      return()
    }
    coverage <- coverages[time_index]
    if(coverage == 0){
      return()
    }
    in_age <- variables$birth$get_index_of(
      timestep - parameters$pmc_ages
    )$to_vector()
    target <- in_age[sample_intervention(in_age, 'pmc', coverage, correlations)]
    
    renderer$render('n_pmc_treated', length(target), timestep)
    treated <- individual::Bitset$new(parameters$human_population)$insert(target)
    
    to_move <- calculate_successful_treatments(
      parameters,
      treated,
      rep(drug, treated$size()),
      timestep,
      renderer,
      "pmc_")
    
    update_mass_drug_admin(
      to_move,
      variables,
      parameters,
      timestep,
      drug
    )
  }
}
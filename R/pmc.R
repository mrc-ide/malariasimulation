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
) {
  function(timestep) {
    timestep_index <- match_timestep(ts = timesteps, t = timestep)
    if(timestep_index == 0){
      return()
    }
    coverage <- coverages[timestep_index]
    if(coverage == 0){
      return()
    }
    
    age <- get_age(variables$birth$get_values(), timestep)
    
    in_age <- which(age %in% parameters$pmc_ages)
    target <- in_age[sample_intervention(in_age, 'pmc', coverage, correlations)]
    
    renderer$render('n_pmc_treated', length(target), timestep)
    
    successful_treatments <- bernoulli(
      length(target),
      parameters$drug_efficacy[[drug]]
    )
    to_move <- individual::Bitset$new(parameters$human_population)
    to_move$insert(target[successful_treatments])
    
    if (to_move$size() > 0) {
      # Move Diseased
      diseased <- variables$state$get_index_of(c('D', 'A'))$and(to_move)
      if (diseased$size() > 0) {
        variables$state$queue_update('Tr', diseased)
      }
      
      # Move everyone else
      other <- to_move$copy()$and(diseased$not(TRUE))
      if (other$size() > 0) {
        variables$state$queue_update('S', other)
      }
      
      # Update infectivity
      variables$infectivity$queue_update(
        variables$infectivity$get_values(
          to_move
        ) * parameters$drug_rel_c[[drug]],
        to_move
      )
      
      # Update drug
      variables$drug$queue_update(drug, to_move)
      variables$drug_time$queue_update(timestep, to_move)
    }
  }
}
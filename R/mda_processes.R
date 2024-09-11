#' @title Create listeners for MDA events
#' @param variables the variables available in the model
#' @param drug the drug to administer
#' @param timesteps timesteps for each round
#' @param coverages the coverage for each round
#' @param min_ages minimum age for the target population for each round
#' @param max_ages maximum age for the target population for each round
#' @param correlations correlation parameters
#' @param int_name the name of this intervention (either 'smc' or 'mda')
#' @param parameters the model parameters
#' @param renderer the model renderer object
#' @description will create a listener for administering each round of drugs
#' @noRd
create_mda_listeners <- function(
    variables,
    drug,
    timesteps,
    coverages,
    min_ages,
    max_ages,
    correlations,
    int_name,
    parameters,
    renderer
) {
  
  renderer$set_default(paste0('n_', int_name, '_treated'), 0)
  renderer$set_default(paste0('n_', int_name, '_drug_efficacy_failures'), 0)
  renderer$set_default(paste0('n_', int_name, '_successfully_treated'), 0)
  
  if(parameters$antimalarial_resistance){
    renderer$set_default(paste0('n_', int_name, '_early_treatment_failure'), 0)
    renderer$set_default(paste0('n_', int_name, '_slow_parasite_clearance'), 0)
  }
  
  function(timestep) {
    time_index = which(timesteps == timestep)
    if(time_index == 0){
      return()
    }
    coverage <- coverages[[time_index]]
    if(coverage == 0){
      return()
    }
    in_age <- variables$birth$get_index_of(
      a = timestep - max_ages[[time_index]],
      b = timestep - min_ages[[time_index]]
    )$to_vector()
    target <- in_age[sample_intervention(in_age, int_name, coverage, correlations)]
    
    renderer$render(paste0('n_', int_name, '_treated'), length(target), timestep)
    treated <- individual::Bitset$new(parameters$human_population)$insert(target)
    
    to_move <- calculate_successful_treatments(
      parameters,
      treated,
      rep(drug, treated$size()),
      timestep,
      renderer,
      paste0(int_name,"_")
      )
    
    update_mass_drug_admin(
      to_move,
      variables,
      parameters,
      timestep,
      drug
    )
    
  }
}

#' @title Update individuals during MDA/PMC
#' @description Updates individuals disease states, infectivity, dt and drug variables
#' @param target a list containing the successfully treated, the drug used and resistance parameters
#' @param variables the variables available in the model
#' @param parameters the model parameters
#' @param timestep the current timestep
#' @param drug the drug to administer
#' @noRd
update_mass_drug_admin <- function(
    target,
    variables,
    parameters,
    timestep,
    drug
){
  
  if (target$successfully_treated$size() > 0) {
    # Move clinical and detectable asymptomatic into treated
    clinical <- variables$state$get_index_of('D')
    asymptomatic <- variables$state$get_index_of('A')
    detectable <- calculate_asymptomatic_detectable(
      variables$state,
      variables$birth,
      variables$id,
      parameters,
      timestep
    )
    to_treated <- clinical$or(asymptomatic$and(detectable))$and(target$successfully_treated)
    
    if(parameters$antimalarial_resistance) {
      dt_update_vector <- target$dt_spc_combined[
        target$successfully_treated$copy()$and(to_treated)$to_vector()
      ]
    } else {
      dt_update_vector <- parameters$dt
    }
    
    update_infection(
      variables$state,
      'Tr',
      variables$infectivity,
      variables$infectivity$get_values(to_treated) * parameters$drug_rel_c[[drug]],
      variables$progression_rates,
      1/dt_update_vector,
      to_treated
    )
    
    # Move everyone else (susceptible, subpatent, non-detected asymptomatic and treated) to susceptible
    other <- target$successfully_treated$copy()$and(to_treated$not(TRUE))
    if (other$size() > 0) {
      update_infection(
        variables$state,
        "S",
        variables$infectivity,
        0,
        variables$progression_rates,
        0,
        other
      )
    }
    
    # Update drug
    variables$drug$queue_update(drug, target$successfully_treated)
    variables$drug_time$queue_update(timestep, target$successfully_treated)
    
  }
}

#' @title Calculate asymptomatic detectable individuals
#'
#' @description Sample a bitset of individuals who are asymptomatic and also
#' detectable by microscopy
#' @param state human infection state
#' @param birth variable for birth of the individual
#' @param immunity to detection
#' @param parameters model parameters
#' @param timestep current timestep
#'
#' @noRd
calculate_asymptomatic_detectable <- function(
    state,
    birth,
    immunity,
    parameters,
    timestep
) {
  asymptomatic <- state$get_index_of('A')
  prob <- probability_of_detection(
    get_age(birth$get_values(asymptomatic), timestep),
    immunity$get_values(asymptomatic),
    parameters
  )
  bitset_at(asymptomatic, bernoulli_multi_p(prob))
}

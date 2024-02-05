#' @title Parameterise antimalarial resistance
#'
#' @param parameters the model parameters
#' @param drug the index of the drug which resistance is being set, as set by the set_drugs() function, in the parameter list
#' @param timesteps vector of time steps for each update to resistance proportion and resistance outcome probability
#' @param artemisinin_resistance vector of updates to the proportions of infections that are artemisinin resistant at time t
#' @param partner_drug_resistance vector of updates to the proportions of infections that are partner-drug resistant at time t
#' @param slow_parasite_clearance_prob vector of updates to the proportion of artemisinin-resistant infections that result in early treatment failure
#' @param early_treatment_failure_prob vector of updates to the proportion of artemisinin-resistant infections that result in slow parasite clearance
#' @param late_clinical_failure_prob vector of updates to the proportion of partner-drug-resistant infections that result in late clinical failure
#' @param late_parasitological_prob vector of updates to the proportion of partner-drug-resistant infections that result in late parasitological failure
#' @param reinfection_prob vector of updates to the proportion of partner-drug-resistant infections that result in reinfection during prophylaxis
#' @param slow_parasite_clearance_time single value representing the mean time individual's experiencing slow parasite clearance reside in the treated state
#' @export
set_antimalarial_resistance <- function(parameters,
                                        drug,
                                        timesteps,
                                        artemisinin_resistance,
                                        partner_drug_resistance,
                                        slow_parasite_clearance_prob,
                                        early_treatment_failure_prob,
                                        late_clinical_failure_prob,
                                        late_parasitological_prob,
                                        reinfection_prob,
                                        slow_parasite_clearance_time) {
  
  # Check that the number of values input is equal for each resistance parameter
  if(any(c(length(artemisinin_resistance),
           length(partner_drug_resistance), 
           length(slow_parasite_clearance_prob), 
           length(early_treatment_failure_prob), 
           length(late_clinical_failure_prob), 
           length(late_parasitological_prob), 
           length(reinfection_prob)) != length(timesteps))) {
    stop("Length of one or more resistance parameter vectors does not match time steps specified for update")
  }
  
  # Ensure resistance proportions bounded between 0 and 1:
  if(any(artemisinin_resistance < 0 | artemisinin_resistance > 1 |
         partner_drug_resistance < 0 | partner_drug_resistance > 1)) {
    stop("Artemisinin and partner-drug resistance proportions must fall between 0 and 1")
  }
  
  # Ensure resistance outcome probabilities bounded between 0 and 1:
  if(any(slow_parasite_clearance_prob < 0 | slow_parasite_clearance_prob > 1 |
         early_treatment_failure_prob < 0 | early_treatment_failure_prob > 1 |
         late_clinical_failure_prob < 0 | late_clinical_failure_prob > 1 |
         late_parasitological_prob < 0 | late_parasitological_prob > 1 |
         reinfection_prob < 0 | reinfection_prob > 1)) {
    stop("Resistance outcome probabilities must fall between 0 and 1")
  }
  
  # Ensure slow_parasite_clearance_time of length 1 and is positive:
  if(slow_parasite_clearance_time <= 0) {
    stop("Error: slow_parasite_clearance_time is non-positive")
  }
  
  # Set antimalarial_resistance to TRUE:
  parameters$antimalarial_resistance <- TRUE
  
  # Store the number of drugs for which parameters are available in the parameter list:
  n_drugs <- length(parameters$drug_efficacy)
  
  # If the drug index falls outside range 1:n_drugs, terminate the operation:
  if (drug < 1 | drug > n_drugs) {
    stop('Drug index is invalid, please set drugs using set_drugs')
  }
  
  # Check the drug_index for the drug we're setting parameters for:
  drug_index <- which(parameters$antimalarial_resistance_drug == drug)
  
  # If drug_index not already assigned, assign the drug the next available index:
  if (length(drug_index) == 0) {
    drug_index <- length(parameters$antimalarial_resistance_drug) + 1
  }
  
  # Append the drug for which resistance is being assigned to the generated index:
  parameters$antimalarial_resistance_drug[[drug_index]] <- drug
  
  # Append the timesteps on which the resistance proportions are to be updated:
  parameters$antimalarial_resistance_timesteps[[drug_index]] <- timesteps
  
  # Append the proportions of all malarial infections that are artemisinin or partner-drug
  # resistant:
  parameters$prop_artemisinin_resistant[[drug_index]] <- artemisinin_resistance
  parameters$prop_partner_drug_resistant[[drug_index]] <- partner_drug_resistance
  
  # Append the probabilities that individuals will experience the resistance outcomes:
  parameters$slow_parasite_clearance_prob[[drug_index]] <- slow_parasite_clearance_prob
  parameters$early_treatment_failure_prob[[drug_index]] <- early_treatment_failure_prob
  parameters$late_clinical_failure_prob[[drug_index]] <- late_clinical_failure_prob
  parameters$late_parasitological_failure_prob[[drug_index]] <- late_parasitological_prob
  parameters$reinfection_during_prophylaxis[[drug_index]] <- reinfection_prob
  
  # Append the slow parasite clearance time:
  parameters$dt_slow_parasite_clearance[[drug_index]] <- slow_parasite_clearance_time
  
  # Return the parameter list:
  parameters
  
}

#' @title Retrieve resistance parameters
#' @description
#' Retrieve the resistance parameters for each individual that has received clinical treatment given
#' the drug they have been administered and the resistance status in the current time step
#' 
#' @param parameters the model parameters
#' @param drug vector of integers representing the drugs administered to each individual seeking treatment
#' @param timestep the current time step
get_antimalarial_resistance_parameters <- function(parameters, drugs, timestep) {
  
  # Check that resistance has been switched on:
  if(parameters$antimalarial_resistance != TRUE) {
    stop("Error: Antimalarial resistance has not been parameterised; antimalarial_resistance = FALSE")
  }
  
  # Retrieve the resistance parameter index of the drug administered for each individual:
  antimalarial_resistance_drug_index <- as.numeric(parameters$antimalarial_resistance_drug)[drugs]
  
  # Establish vectors to store the artemisinin and resistance proportions associated with the drug
  # administered to each individual:
  artemisinin_resistance_proportion <- vector()
  partner_drug_resistance_proportion <- vector()
  
  # Establish the vectors to store the resistance outcome probabilities associated with the drug
  # administered to each individual:
  slow_parasite_clearance_probability <- vector()
  early_treatment_failure_probability <- vector()
  late_clinical_failure_probability <- vector()
  late_parasitological_failure_probability <- vector()
  reinfection_during_prophylaxis_probability <- vector()
  
  # Establish a vector to store the infectious state durations associated with resistance:
  dt_slow_parasite_clearance <- vector()
  
  # Loop through and assign assign the correct antimalarial resistance parameters for each individual
  # given the drug they were administered:
  for(i in seq_along(parameters$antimalarial_resistance_drug)) {
    
    # Retrieve the indices of individuals administered drug i
    drug_indices <- which(drugs == i)
    
    # Determine the index of drug i in the antimalarial resistance parameters
    drug_index <- which(parameters$antimalarial_resistance_drug == i)
    
    # Determine which antimalarial resistance parameters apply in the current time step for drug i:
    resistance_timestep <- match_timestep(ts = parameters$antimalarial_resistance_timesteps[[drug_index]], t = timestep)
    
    # Retrieve and assign the resistance proportions for drug i:
    artemisinin_resistance_proportion[drug_indices] <- parameters$prop_artemisinin_resistant[[drug_index]][resistance_timestep]
    partner_drug_resistance_proportion[drug_indices] <- parameters$prop_partner_drug_resistant[[drug_index]][resistance_timestep]
    
    # Retrieve and assign the resistance outcome probabilities for drug i:
    slow_parasite_clearance_probability[drug_indices] <- parameters$slow_parasite_clearance_prob[[drug_index]][resistance_timestep]
    early_treatment_failure_probability[drug_indices] <- parameters$early_treatment_failure_prob[[drug_index]][resistance_timestep]
    late_clinical_failure_probability[drug_indices] <- parameters$late_clinical_failure_prob[[drug_index]][resistance_timestep]
    late_parasitological_failure_probability[drug_indices] <- parameters$late_parasitological_failure_prob[[drug_index]][resistance_timestep]
    reinfection_during_prophylaxis_probability[drug_indices] <- parameters$reinfection_during_prophylaxis[[drug_index]][resistance_timestep]
    
    # Retrieve and assign the infectious state durations associated with resistance to drug i:
    dt_slow_parasite_clearance[drug_indices] <- parameters$dt_slow_parasite_clearance[[drug_index]]
  }
  
  # Store the vectors of parameters in a list:
  resistance_parameters <- list()
  resistance_parameters$artemisinin_resistance_proportion <- artemisinin_resistance_proportion
  resistance_parameters$partner_drug_resistance_proportion <- partner_drug_resistance_proportion
  resistance_parameters$slow_parasite_clearance_probability <- slow_parasite_clearance_probability
  resistance_parameters$early_treatment_failure_probability <- early_treatment_failure_probability
  resistance_parameters$late_clinical_failure_probability <- late_clinical_failure_probability
  resistance_parameters$late_parasitological_failure_probability <- late_parasitological_failure_probability
  resistance_parameters$reinfection_during_prophylaxis_probability <- reinfection_during_prophylaxis_probability
  resistance_parameters$dt_slow_parasite_clearance <- dt_slow_parasite_clearance
  
  # Check that each resistance parameter entry contains the correct number of values:
  if(any(c(length(resistance_parameters$artemisinin_resistance_proportion),
           length(resistance_parameters$partner_drug_resistance_proportion), 
           length(resistance_parameters$slow_parasite_clearance_probability), 
           length(resistance_parameters$early_treatment_failure_probability), 
           length(resistance_parameters$late_clinical_failure_probability), 
           length(resistance_parameters$late_parasitological_failure_probability), 
           length(resistance_parameters$reinfection_during_prophylaxis_probability),
           length(resistance_parameters$dt_slow_parasite_clearance)) != length(drugs))) {
    stop("Length of one or more resistance parameter vectors does not match length of drugs vector")
  }
  
  # Return the list of parameters:
  return(resistance_parameters)
  
}

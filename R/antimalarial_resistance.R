#' @title Parameterise antimalarial resistance
#' @description
#' Parameterise antimalarial resistance
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
  
  if(any(c(length(artemisinin_resistance),
           length(partner_drug_resistance), 
           length(slow_parasite_clearance_prob), 
           length(early_treatment_failure_prob), 
           length(late_clinical_failure_prob), 
           length(late_parasitological_prob), 
           length(reinfection_prob)) != length(timesteps))) {
    stop("Length of one or more resistance parameter vectors does not match time steps specified for update")
  }
  
  if(any(artemisinin_resistance < 0 | artemisinin_resistance > 1 |
         partner_drug_resistance < 0 | partner_drug_resistance > 1)) {
    stop("Artemisinin and partner-drug resistance proportions must fall between 0 and 1")
  }
  
  if(any(slow_parasite_clearance_prob < 0 | slow_parasite_clearance_prob > 1 |
         early_treatment_failure_prob < 0 | early_treatment_failure_prob > 1 |
         late_clinical_failure_prob < 0 | late_clinical_failure_prob > 1 |
         late_parasitological_prob < 0 | late_parasitological_prob > 1 |
         reinfection_prob < 0 | reinfection_prob > 1)) {
    stop("Resistance outcome probabilities must fall between 0 and 1")
  }

    if(length(slow_parasite_clearance_time) != 1) {
    stop("Error: length of slow_parasite_clearance_time not equal to 1")
  }
  
  if(slow_parasite_clearance_time <= 0) {
    stop("Error: slow_parasite_clearance_time is non-positive")
  }
  
  parameters$antimalarial_resistance <- TRUE
  
  n_drugs <- length(parameters$drug_efficacy)
  
  if (drug < 1 | drug > n_drugs) {
    stop('Drug index is invalid, please set drugs using set_drugs')
  }
  
  drug_index <- which(parameters$antimalarial_resistance_drug == drug)
  
  if (length(drug_index) == 0) {
    drug_index <- length(parameters$antimalarial_resistance_drug) + 1
  }
  
  parameters$antimalarial_resistance_drug[[drug_index]] <- drug
  parameters$antimalarial_resistance_timesteps[[drug_index]] <- timesteps
  parameters$prop_artemisinin_resistant[[drug_index]] <- artemisinin_resistance
  parameters$prop_partner_drug_resistant[[drug_index]] <- partner_drug_resistance
  parameters$slow_parasite_clearance_prob[[drug_index]] <- slow_parasite_clearance_prob
  parameters$early_treatment_failure_prob[[drug_index]] <- early_treatment_failure_prob
  parameters$late_clinical_failure_prob[[drug_index]] <- late_clinical_failure_prob
  parameters$late_parasitological_failure_prob[[drug_index]] <- late_parasitological_prob
  parameters$reinfection_during_prophylaxis[[drug_index]] <- reinfection_prob
  parameters$dt_slow_parasite_clearance[[drug_index]] <- slow_parasite_clearance_time
  
  return(parameters)
  
}

#' @title Retrieve resistance parameters
#' @description
#' Retrieve the resistance parameters associated with the drug each individual receiving clinical 
#' treatment has been administered in the current time step.
#' 
#' @param parameters the model parameters
#' @param drug vector of integers representing the drugs administered to each individual receiving treatment
#' @param timestep the current time step
get_antimalarial_resistance_parameters <- function(parameters, drugs, timestep) {
  
  if(!parameters$antimalarial_resistance) {
    stop("Error: Antimalarial resistance has not been parameterised; antimalarial_resistance = FALSE")
  }
  
  blank_vector <- numeric(length = length(drugs))
  artemisinin_resistance_proportion <- blank_vector
  partner_drug_resistance_proportion <- blank_vector
  slow_parasite_clearance_probability <- blank_vector
  early_treatment_failure_probability <- blank_vector
  late_clinical_failure_probability <- blank_vector
  late_parasitological_failure_probability <- blank_vector
  reinfection_during_prophylaxis_probability <- blank_vector
  dt_slow_parasite_clearance <- rep(parameters$dt, length = length(drugs))
  
  for(i in seq_along(parameters$antimalarial_resistance_drug)) {
    drug <- parameters$antimalarial_resistance_drug[[i]]
    treated_with_drug <- which(drugs == drug)
    resistance_timestep <- match_timestep(ts = parameters$antimalarial_resistance_timesteps[[i]], t = timestep)
    artemisinin_resistance_proportion[treated_with_drug] <- parameters$prop_artemisinin_resistant[[i]][resistance_timestep]
    partner_drug_resistance_proportion[treated_with_drug] <- parameters$prop_partner_drug_resistant[[i]][resistance_timestep]
    slow_parasite_clearance_probability[treated_with_drug] <- parameters$slow_parasite_clearance_prob[[i]][resistance_timestep]
    early_treatment_failure_probability[treated_with_drug] <- parameters$early_treatment_failure_prob[[i]][resistance_timestep]
    late_clinical_failure_probability[treated_with_drug] <- parameters$late_clinical_failure_prob[[i]][resistance_timestep]
    late_parasitological_failure_probability[treated_with_drug] <- parameters$late_parasitological_failure_prob[[i]][resistance_timestep]
    reinfection_during_prophylaxis_probability[treated_with_drug] <- parameters$reinfection_during_prophylaxis[[i]][resistance_timestep]
    dt_slow_parasite_clearance[treated_with_drug] <- parameters$dt_slow_parasite_clearance[[i]]
  }
  
  resistance_parameters <- list()
  resistance_parameters$artemisinin_resistance_proportion <- artemisinin_resistance_proportion
  resistance_parameters$partner_drug_resistance_proportion <- partner_drug_resistance_proportion
  resistance_parameters$slow_parasite_clearance_probability <- slow_parasite_clearance_probability
  resistance_parameters$early_treatment_failure_probability <- early_treatment_failure_probability
  resistance_parameters$late_clinical_failure_probability <- late_clinical_failure_probability
  resistance_parameters$late_parasitological_failure_probability <- late_parasitological_failure_probability
  resistance_parameters$reinfection_during_prophylaxis_probability <- reinfection_during_prophylaxis_probability
  resistance_parameters$dt_slow_parasite_clearance <- dt_slow_parasite_clearance
  
  return(resistance_parameters)
  
}
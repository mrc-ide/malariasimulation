#' @title Parameterise antimalarial resistance
#' @description
#' Parameterise antimalarial resistance
#' 
#' @param parameters the model parameters
#' @param drug the index of the drug which resistance is being set, as set by the set_drugs() function, in the parameter list
#' @param timesteps vector of time steps for each update to resistance proportion and resistance outcome probability
#' @param artemisinin_resistance_proportion vector of updates to the proportions of infections that are artemisinin resistant at time t
#' @param partner_drug_resistance_proportion vector of updates to the proportions of infections that are partner-drug resistant at time t
#' @param slow_parasite_clearance_probability vector of updates to the proportion of artemisinin-resistant infections that result in early treatment failure
#' @param early_treatment_failure_probability vector of updates to the proportion of artemisinin-resistant infections that result in slow parasite clearance
#' @param late_clinical_failure_probability vector of updates to the proportion of partner-drug-resistant infections that result in late clinical failure
#' @param late_parasitological_failure_probability vector of updates to the proportion of partner-drug-resistant infections that result in late parasitological failure
#' @param reinfection_during_prophylaxis_probability vector of updates to the proportion of partner-drug-resistant infections that result in reinfection during prophylaxis
#' @param slow_parasite_clearance_time single value representing the mean time individual's experiencing slow parasite clearance reside in the treated state
#' @export
set_antimalarial_resistance <- function(parameters,
                                        drug,
                                        timesteps,
                                        artemisinin_resistance_proportion,
                                        partner_drug_resistance_proportion,
                                        slow_parasite_clearance_probability,
                                        early_treatment_failure_probability,
                                        late_clinical_failure_probability,
                                        late_parasitological_failure_probability,
                                        reinfection_during_prophylaxis_probability,
                                        slow_parasite_clearance_time) {
  
  if(any(partner_drug_resistance_proportion > 0,
         late_clinical_failure_probability > 0,
         late_parasitological_failure_probability > 0,
         reinfection_during_prophylaxis_probability > 0)) {
    stop("Parameters set for unimplemented feature - late clinical failure, late parasitological failure, or reinfection during prophylaxis")
  }
  
  if(any(c(length(artemisinin_resistance_proportion),
           length(partner_drug_resistance_proportion), 
           length(slow_parasite_clearance_probability), 
           length(early_treatment_failure_probability), 
           length(late_clinical_failure_probability), 
           length(late_parasitological_failure_probability), 
           length(reinfection_during_prophylaxis_probability)) != length(timesteps))) {
    stop("Length of one or more resistance parameter vectors does not match time steps specified for update")
  }
  
  if(any(artemisinin_resistance_proportion < 0 | artemisinin_resistance_proportion > 1 |
         partner_drug_resistance_proportion < 0 | partner_drug_resistance_proportion > 1)) {
    stop("Artemisinin and partner-drug resistance proportions must fall between 0 and 1")
  }
  
  if(any(slow_parasite_clearance_probability < 0 | slow_parasite_clearance_probability > 1 |
         early_treatment_failure_probability < 0 | early_treatment_failure_probability > 1 |
         late_clinical_failure_probability < 0 | late_clinical_failure_probability > 1 |
         late_parasitological_failure_probability < 0 | late_parasitological_failure_probability > 1 |
         reinfection_during_prophylaxis_probability < 0 | reinfection_during_prophylaxis_probability > 1)) {
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
  parameters$artemisinin_resistance_prop[[drug_index]] <- artemisinin_resistance_proportion
  parameters$partner_drug_resistance_prop[[drug_index]] <- partner_drug_resistance_proportion
  parameters$slow_parasite_clearance_prob[[drug_index]] <- slow_parasite_clearance_probability
  parameters$early_treatment_failure_prob[[drug_index]] <- early_treatment_failure_probability
  parameters$late_clinical_failure_prob[[drug_index]] <- late_clinical_failure_probability
  parameters$late_parasitological_failure_prob[[drug_index]] <- late_parasitological_failure_probability
  parameters$reinfection_during_prophylaxis_probability[[drug_index]] <- reinfection_during_prophylaxis_probability
  parameters$dt_slow_parasite_clearance[[drug_index]] <- slow_parasite_clearance_time
  
  return(parameters)
  
}

#' @title Retrieve resistance parameters
#' @description
#' Retrieve the resistance parameters associated with the drug each individual receiving clinical 
#' treatment has been administered in the current time step.
#' 
#' @param parameters the model parameters
#' @param drugs vector of integers representing the drugs administered to each individual receiving treatment
#' @param timestep the current time step
get_antimalarial_resistance_parameters <- function(parameters, drugs, timestep) {
  
  if(!parameters$antimalarial_resistance) {
    stop("Error: Antimalarial resistance has not been parameterised; antimalarial_resistance = FALSE")
  }
  
  blank_vector <- numeric(length = length(drugs))
  artemisinin_resistance_prop <- blank_vector
  partner_drug_resistance_prop <- blank_vector
  slow_parasite_clearance_prob <- blank_vector
  early_treatment_failure_prob <- blank_vector
  late_clinical_failure_prob <- blank_vector
  late_parasitological_failure_prob <- blank_vector
  reinfection_during_prophylaxis_prob <- blank_vector
  dt_slow_parasite_clearance <- rep(parameters$dt, length = length(drugs))
  
  for(i in seq_along(parameters$antimalarial_resistance_drug)) {
    drug <- parameters$antimalarial_resistance_drug[[i]]
    treated_with_drug <- which(drugs == drug)
    resistance_timestep <- match_timestep(ts = parameters$antimalarial_resistance_timesteps[[i]], t = timestep)
    artemisinin_resistance_prop[treated_with_drug] <- parameters$artemisinin_resistance_prop[[i]][resistance_timestep]
    partner_drug_resistance_prop[treated_with_drug] <- parameters$partner_drug_resistance_prop[[i]][resistance_timestep]
    slow_parasite_clearance_prob[treated_with_drug] <- parameters$slow_parasite_clearance_prob[[i]][resistance_timestep]
    early_treatment_failure_prob[treated_with_drug] <- parameters$early_treatment_failure_prob[[i]][resistance_timestep]
    late_clinical_failure_prob[treated_with_drug] <- parameters$late_clinical_failure_prob[[i]][resistance_timestep]
    late_parasitological_failure_prob[treated_with_drug] <- parameters$late_parasitological_failure_prob[[i]][resistance_timestep]
    reinfection_during_prophylaxis_prob[treated_with_drug] <- parameters$reinfection_during_prophylaxis_prob[[i]][resistance_timestep]
    dt_slow_parasite_clearance[treated_with_drug] <- parameters$dt_slow_parasite_clearance[[i]]
  }
  
  resistance_parameters <- list()
  resistance_parameters$artemisinin_resistance_prop <- artemisinin_resistance_prop
  resistance_parameters$partner_drug_resistance_prop <- partner_drug_resistance_prop
  resistance_parameters$slow_parasite_clearance_prob <- slow_parasite_clearance_prob
  resistance_parameters$early_treatment_failure_prob <- early_treatment_failure_prob
  resistance_parameters$late_clinical_failure_prob <- late_clinical_failure_prob
  resistance_parameters$late_parasitological_failure_prob <- late_parasitological_failure_prob
  resistance_parameters$reinfection_during_prophylaxis_prob <- reinfection_during_prophylaxis_prob
  resistance_parameters$dt_slow_parasite_clearance <- dt_slow_parasite_clearance
  
  return(resistance_parameters)
  
}

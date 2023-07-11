#' @title Parameterise antimalarial resistance
#'
#' @param parameters the model parameters
#' @param drug the index of the drug which resistance is being set for in the parameter list
#' @param timesteps vector of timesteps for each update to resistance proportion and/or phenotype probability
#' @param artemisinin_resistance vector of updates to the proportions of infections that are artemisinin resistant at time t
#' @param partner_drug_resistance vector of updates to the proportions of infections that are partner-drug resistant at time t
#' @param slow_parasite_clearance_prob vector of updates to the proportion of artemisinin-resistant infections that result in early treatment failure
#' @param early_treatment_failure_prob vector of updates to the proportion of artemisinin-resistant infections that result in slow parasite clearance
#' @param late_clinical_failure_prob vector of updates to the proportion of partner-drug-resistant infections that result in late clinical failure
#' @param late_parasitological_prob vector of updates to the proportion of partner-drug-resistant infections that result in late parasitologica; failure
#' @param reinfection_prob vector of updates to the proportion of partner-drug-resistant infections that result in reinfection during prophylaxis
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
                                        reinfection_prob) {
  
  # Check that the number of values input is equal for each resistance parameter
  if(length(artemisinin_resistance) != length(timesteps) |
     length(partner_drug_resistance) != length(timesteps) |
     length(slow_parasite_clearance_prob) != length(timesteps) |
     length(early_treatment_failure_prob) != length(timesteps) |
     length(late_clinical_failure_prob) != length(timesteps) |
     length(late_parasitological_prob) != length(timesteps) |
     length(reinfection_prob) != length(timesteps)) {
    stop("Number of resistance parameter vectors do not match time steps specified for update")
  } else {
    print("OK")
  }
  
  # Check that the proportion of people with artemisinin and partner-drug resistance
  # are bounded between 0 and 1: 
  for(i in length(artemisinin_resistance)) {
    if(artemisinin_resistance[i] < 0 | artemisinin_resistance[i] > 1 |
       partner_drug_resistance[i] < 0 | partner_drug_resistance[i] > 1) {
      stop("Artemisinin and partner-drug resistance proportions must fall between 0 and 1")
    } else {
      print("OK")
    }
  }
  
  # Ensure resistance phenotype probabilities bounded between 0 and 1:
  for(i in 1:length(slow_parasite_clearance_prob)) {
    if(slow_parasite_clearance_prob[i] < 0 | slow_parasite_clearance_prob[i] > 1 |
       early_treatment_failure_prob[i] < 0 | early_treatment_failure_prob[i] > 1 |
       late_clinical_failure_prob[i] < 0 | late_clinical_failure_prob[i] > 1 |
       late_parasitological_prob[i] < 0 | late_parasitological_prob[i] > 1 |
       reinfection_prob[i] < 0 | reinfection_prob[i] > 1) {
      stop("Resistance phenotype probabilities must fall between 0 and 1")
    } else {
      print("OK")
    }
  }
  
  # Set antimalarial_resistance to TRUE
  parameters$antimalarial_resistance <- TRUE
  
  # Store the number of drugs for which parameters are available in the parameter list:
  n_drugs <- length(parameters$drug_efficacy)
  
  # If the drug index falls outside range 1:n_drugs, terminate the operation:
  if (drug < 1 | drug > n_drugs) {
    stop('Drug index is invalid, please set drugs using set_drugs')
  } else {
    print("OK")
  }
  
  # Check the drug_index for the drug setting parameters for
  drug_index <- which(parameters$antimalarial_resistance_drug == drug)
  
  # If drug_index is currently unpopulated
  if (length(drug_index) == 0) {
    drug_index <- length(parameters$antimalarial_resistance_drug) + 1
  }
  
  # Append the drug for which resistance is being assigned
  parameters$antimalarial_resistance_drug[[drug_index]] <- drug
  
  # Append the timesteps on which the resistance proportions are to be updated
  parameters$antimalarial_resistance_timesteps[[drug_index]] <- timesteps
  
  # Append the proportions of all malarial infections that are artemisinin or partner-drug
  # resistant
  parameters$prop_artemisinin_resistant[[drug_index]] <- artemisinin_resistance
  parameters$prop_partner_drug_resistant[[drug_index]] <- partner_drug_resistance
  
  # Append the probabilities that
  parameters$slow_parasite_clearance_prob[[drug_index]] <- slow_parasite_clearance_prob
  parameters$early_treatment_failure_prob[[drug_index]] <- early_treatment_failure_prob
  parameters$late_clinical_failure_prob[[drug_index]] <- late_clinical_failure_prob
  parameters$late_parasitological_failure_prob[[drug_index]] <- late_parasitological_prob
  parameters$reinfection_during_prophylaxis[[drug_index]] <- reinfection_prob
  
  # Return the parameter list:
  parameters
  
}
#' @title Preset parameters for the DHA-PQP drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a vector of preset parameters for DHA-PQP drug
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.09434, drug_prophylaxis_shape: 4.4, drug_prophylaxis_scale: 28.1
#' @export
DHA_PQP_params <- c(.95, 0.09434, 4.4, 28.1)

#' @title Preset parameters for the AL drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @export
AL_params <- c(.95, 0.05094, 11.3, 10.6)

#' @title Preset parameters for the SP-AQ drug
#' @export
SP_AQ_params <- c(0.9, 0.32, 4.3, 38.1)

#' @title Parameterise drugs to use in the model
#'
#' @param parameters the model parameters
#' @param drugs a list of drug parameters, can be set using the above presets
#' @export
set_drugs <- function(parameters, drugs) {
  keys <- c(
    'drug_efficacy',
    'drug_rel_c',
    'drug_prophylaxis_shape',
    'drug_prophylaxis_scale'
  )
  for (drug in seq_along(drugs)) {
    for (i in seq_along(drugs[[drug]])) {
      parameters[[keys[[i]]]] <- c(parameters[[keys[[i]]]], drugs[[drug]][[i]])
    }
  }
  parameters
}

#' @title Parameterise clinical treatment
#'
#' @param parameters the model parameters
#' @param drug the index of the drug to set as a treatment
#' @param timesteps vector of timesteps for each coverage change
#' @param coverages vector of coverages for this drug
#' @export
set_clinical_treatment <- function(parameters, drug, timesteps, coverages) {
  n_drugs <- length(parameters$drug_efficacy)
  if (drug < 1 | drug > n_drugs) {
    stop('Drug index is invalid, please set drugs using set_drugs')
  }
  drug_index <- which(parameters$clinical_treatment_drugs == drug)
  if (length(drug_index) == 0) {
    drug_index <- length(parameters$clinical_treatment_drugs) + 1
  }
  parameters$clinical_treatment_drugs[[drug_index]] <- drug
  parameters$clinical_treatment_timesteps[[drug_index]] <- timesteps
  parameters$clinical_treatment_coverages[[drug_index]] <- coverages
  last_timestep <- max(unlist(parameters$clinical_treatment_timesteps))

  for (t in seq(last_timestep)) {
    if (sum(get_treatment_coverages(parameters, t)) > 1) {
      stop('The sum of drug coverages cannot be greater than 1 at any timestep')
    }
  }
  parameters
}

get_treatment_coverages <- function(parameters, timestep) {
  vnapply(
    seq_along(parameters$clinical_treatment_drugs),
    function(drug_index) {
      previous <- which(
        parameters$clinical_treatment_timesteps[[drug_index]] <= timestep
      )
      if (length(previous) == 0) {
        return(0)
      }
      last_set <- max(previous)
      parameters$clinical_treatment_coverages[[drug_index]][[last_set]]
    }
  )
}

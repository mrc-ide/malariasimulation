#' @title Preset parameters for the DHA-PQP drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a vector of preset parameters for the DHA-PQP drug (dihydroartemisinin-piperaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.09434, drug_prophylaxis_shape: 4.4, drug_prophylaxis_scale: 28.1
#' @export
DHA_PQP_params <- c(.95, 0.09434, 4.4, 28.1, "falciparum")

#' @title Preset parameters for the AL drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a vector of preset parameters for the AL drug (artemether-lumefantrine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.05094, drug_prophylaxis_shape: 11.3, drug_prophylaxis_scale: 10.6
#' @export
AL_params <- c(.95, 0.05094, 11.3, 10.6, "falciparum")

#' @title Preset parameters for the SP-AQ drug
#' @details Use a vector of preset parameters for the SP-AQ drug (sulphadoxine-pyrimethamine and amodiaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.9, drug_rel_c: 0.32, drug_prophylaxis_shape: 4.3, drug_prophylaxis_scale: 38.1
#' @export
SP_AQ_params <- c(0.9, 0.32, 4.3, 38.1, "falciparum")

#' @title Preset parameters for the CQ drug
#' @description From SI of Nekkab et al., DOI: 10.1371/journal.pmed.1003535 (2021)
#' @details Use a vector of preset parameters for the CQ drug (chloroquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.899, drug_rel_c: 1, drug_prophylaxis_shape: 1, drug_prophylaxis_scale: 1/28
#' @export
CQ_params <- c(0.899, 1, 1, 1/28, "vivax")

#' @title Preset parameters for the CQ-PQ drug
#' @description From SI of Nekkab et al., DOI: 10.1371/journal.pmed.1003535 (2021)
#' @details Use a vector of preset parameters for the CQ-PQ drug (chloroquine and primaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.948, drug_rel_c: 1, drug_prophylaxis_shape: 1, drug_prophylaxis_scale: 1/28, drug_hypnozoite_efficacy: 0.713, drug_hypnozoite_prophylaxis: 8 days
#' @export
CQ_PQ_params <- c(0.948, 1, 1, 1/28, "vivax", 0.713, 8)

#' @title Preset parameters for the CQ-TQ drug
#' @description From SI of Nekkab et al., DOI: 10.1038/s41467-018-05860-8 (2018)
#' @details Use a vector of preset parameters for the CQ-TQ drug (chloroquine and tafenoquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 1.00, drug_rel_c: 1, drug_prophylaxis_shape: 1, drug_prophylaxis_scale: 1/45, drug_hypnozoite_efficacy: 0.713, drug_hypnozoite_prophylaxis: 45 days
#' @export
CQ_TQ_params <- c(1, 1, 1, 1/45, "vivax", 0.713, 45)

#' @title Preset parameters for the White et a., ACT drug
#' @description From SI of White et al., DOI: 10.1371/journal.pmed.1003535 (2021)
#' @details Use a vector of preset parameters for the ACT drug (Artemisinin combination therapy)
#' @details Default parameters, from L to R, are: drug_efficacy: 1, drug_rel_c: 1, drug_prophylaxis_shape: 1, drug_prophylaxis_scale: 1/14
#' @export
white_ACT_params <- c(0.899, 1, 1, 1/14, "vivax")
# AL_params_vivax <- c(0.899, 1, 1, 1/14, "vivax")

#' @title Preset parameters for the ACT-PQ drug
#' @description From SI of White et al., DOI: 10.1371/journal.pmed.1003535 (2021)
#' @details Use a vector of preset parameters for the ACT-PQ drug (Artemisinin combination therapy and primaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 1, drug_rel_c: 1, drug_prophylaxis_shape: 1, drug_prophylaxis_scale: 1/14, drug_hypnozoite_efficacy: 0.7, drug_hypnozoite_prophylaxis: 14 days
#' @export
white_ACT_PQ_params <- c(1, 1, 1, 1/14, "vivax", 0.7, 14)

#' @title Preset parameters for the ChlQ drug
#' @description From git repo DOI: 10.1371/journal.pmed.1003535 (2021)
#' @details Use a vector of preset parameters for the ChlQ drug (chloroquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 1, drug_rel_c: 1, drug_prophylaxis_shape: 1, drug_prophylaxis_scale: 1/60, drug_hypnozoite_efficacy: 1, drug_hypnozoite_prophylaxis: 60 days
#' @export
white_ACT_TQ_params <- c(1, 1, 1, 1/60, "vivax", 1, 60)


#' @title Parameterise drugs to use in the model
#'
#' @param parameters the model parameters
#' @param drugs a list of drug parameters, can be set using presets
#' @export
set_drugs <- function(parameters, drugs) {
  stopifnot(is.list(drugs))
  keys <- c(
    'drug_efficacy',
    'drug_rel_c',
    'drug_prophylaxis_shape',
    'drug_prophylaxis_scale',
    'parasite_target',
    'drug_hypnozoite_efficacy',
    'drug_hypnozoite_prophylaxis'
  )
  for (drug in seq_along(drugs)) {
    stopifnot("parasite and drug parameters do not align" = parameters$parasite == drugs[[drug]][5])
    for (i in seq_along(drugs[[drug]])[-5]) {
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
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
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

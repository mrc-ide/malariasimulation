#' @title Preset parameters for the DHA-PQP drug (P. falciparum)
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a vector of preset parameters for the DHA-PQP drug (dihydroartemisinin-piperaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.09434, drug_prophylaxis_shape: 4.4, drug_prophylaxis_scale: 28.1
#' @export
DHA_PQP_params <- c(.95, 0.09434, 4.4, 28.1)

#' @title Preset parameters for the AL drug (P. falciparum)
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @details Use a vector of preset parameters for the AL drug (artemether-lumefantrine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.95, drug_rel_c: 0.05094, drug_prophylaxis_shape: 11.3, drug_prophylaxis_scale: 10.6
#' @export
AL_params <- c(.95, 0.05094, 11.3, 10.6)

#' @title Preset parameters for the SP-AQ drug (P. falciparum)
#' @details Use a vector of preset parameters for the SP-AQ drug (sulphadoxine-pyrimethamine and amodiaquine)
#' @details Default parameters, from L to R, are: drug_efficacy: 0.9, drug_rel_c: 0.32, drug_prophylaxis_shape: 4.3, drug_prophylaxis_scale: 38.1
#' @export
SP_AQ_params <- c(0.9, 0.32, 4.3, 38.1)

#' @title Preset parameters for the CQ drug (P. vivax)
#' @description Efficacy from SI of Nekkab et al., DOI: 10.1371/journal.pmed.1003535 (2021),
#' shape and scale consistent with a longer prophylaxis time: 28 days, decreasing gradually
#' @details Use a vector of preset parameters for the CQ drug (chloroquine) acting on P. vivax
#' @details Default parameters, from L to R, are: drug_efficacy: 0.899, drug_rel_c: 0.5, drug_prophylaxis_shape: 20, drug_prophylaxis_scale: 5
#' @export
CQ_params_vivax <- c(0.899, 0.5, 5, 20, 0, 0, 0)


#' @title Preset parameters for the CQ-PQ drug (P. vivax)
#' @description Efficacy from SI of Nekkab et al., DOI: 10.1371/journal.pmed.1003535 (2021),
#' BS shape and scale consistent with a longer prophylaxis time: 28 days, decreasing gradually,
#' LS shape are scale consistent with a 7 day PQ treatment, with rapid decay.
#' @details Use a vector of preset parameters for the CQ-PQ drug (chloroquine and primaquine) acting on P. vivax
#' @details Default parameters, from L to R, are: drug_efficacy: 0.948, drug_rel_c: 0.5, drug_prophylaxis_shape: 5, drug_prophylaxis_scale: 20,
#' drug_hypnozoite_efficacy: 0.713, drug_hypnozoite_prophylaxis_shape: 10, drug_hypnozoite_prophylaxis_scale: 5.5
#' @export
CQ_PQ_params_vivax <- c(0.948, 0.5, 5, 20, 0.713, 10, 5.5)

#' @title Preset parameters for the CQ-TQ drug (P. vivax)
#' @description Efficacy from SI of Nekkab et al., DOI: 10.1371/journal.pmed.1003535 (2021),
#' BS shape and scale consistent with a longer prophylaxis time: 28 days, decreasing gradually,
#' LS shape are scale consistent with a single TQ treatment, with longer prophylaxis: 45 days.
#' @details Use a vector of preset parameters for the CQ-TQ drug (chloroquine and tafenoquine) acting on P. vivax
#' @details Default parameters, from L to R, are: drug_efficacy: 1, drug_rel_c: 0.5, drug_prophylaxis_shape: 5, drug_prophylaxis_scale: 20,
#' drug_hypnozoite_efficacy: 0.713, drug_hypnozoite_prophylaxis_shape: 5.5, drug_hypnozoite_prophylaxis_scale: 30
#' @export
CQ_TQ_params_vivax <- c(1, 0.5, 5, 20, 0.713, 5, 30)

#' @title Parameterise drugs to use in the model
#'
#' @param parameters the model parameters
#' @param drugs a list of drug parameters, can be set using presets
#' @export
set_drugs <- function(parameters, drugs) {
  stopifnot("is.list(drugs) is not TRUE" = is.list(drugs))
  keys <- c(
    'drug_efficacy',
    'drug_rel_c',
    'drug_prophylaxis_shape',
    'drug_prophylaxis_scale',
    # hypnozoite parameters
    'drug_hypnozoite_efficacy',
    'drug_hypnozoite_prophylaxis_shape',
    'drug_hypnozoite_prophylaxis_scale'
  )
  
  if(parameters$parasite == "falciparum"){
    for (drug in seq_along(drugs)) {
      if(length(drugs[[drug]]) != 4){
        warning(paste0("Drug ", drug, " has incorrect number of P. falciparum drug parameters. The number of parameters should be 4."),
                call. = FALSE)
      }
    }
  } else if (parameters$parasite == "vivax"){
    for (drug in seq_along(drugs)) {
      if(length(drugs[[drug]]) != 7){
        warning(paste0("Drug ", drug, " has incorrect number of P. vivax drug parameters. The number of parameters should be 7 for radical cure. To assign a blood stage drug only, set the liver stage drug parameters to 0: see CQ_params_vivax for an example."),
                call. = FALSE)
      }
    }
  }
  
  for (drug in seq_along(drugs)) {
    for (i in seq_along(drugs[[drug]])) {
      parameters[[keys[[i]]]][drug] <- drugs[[drug]][[i]]
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

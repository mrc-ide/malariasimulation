#' @title Preset parameters for the DHC-PQP drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @export
DHC_PQP_params <- c(.95, 0.09434, 4.4, 28.1)

#' @title Preset parameters for the AL drug
#' @description From SI of Commun. 5:5606 doi: 10.1038/ncomms6606 (2014)
#' @export
AL_params <- c(.95, 0.05094, 11.3, 10.6)

#' @title Parameterise drugs to use in the model
#'
#' @param parameters the model parameters
#' @param drugs a list of drug parameters, either preset or created by `create_drug`
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
#' @param ft probability of seeking treatment
#' @param drugs vector of drugs indecies that are available
#' @param coverages vector of coverages for those drugs
#' @export
set_clinical_treatment <- function(parameters, ft, drugs, coverages) {
  parameters$ft <- ft
  if (any(drugs < 1 | drugs > length(parameters$drug_efficacy))) {
    stop('Drug indices are invalid')
  }
  parameters$clinical_treatment_drugs <- drugs
  if (!approx_sum(coverages, 1)) {
    stop('Drug coverages do not sum to 1')
  }
  parameters$clinical_treatment_coverages <- coverages
  parameters
}

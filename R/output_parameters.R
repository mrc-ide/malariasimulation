#' @title Parameterise age grouped output rendering
#'
#' @details this function produces contiguous age groups, inclusive of the lower
#' age limit and exclusive of the upper age limit: e.g., c(0, 10, 100) will produce
#' two age groups: 0-9 and 10-99 in days.
#' @param parameters the model parameters
#' @param age_group age breaks for population size outputs; default = NULL
#' @param incidence age breaks for incidence outputs (pf: D+Tr+A, pv: D+Tr+A+U); default = NULL
#' @param patent_incidence age breaks for patent incidence outputs (LM detectable), (p.v only); default = NULL
#' @param clinical_incidence age breaks for clinical incidence outputs (symptomatic); default = c(0, 1825)
#' @param severe_incidence age breaks for severe incidence outputs (p.f only); default = NULL
#' @param prevalence age breaks for clinical prevalence outputs (pcr and lm detectable infections); default = c(730, 3650)
#' @param hypnozoite_prevalence age breaks for hypnozoite prevalence outputs (p.v only); default = NULL
#' @export
#'
set_epi_outputs <- function(parameters,
                            age_group = NULL,
                            incidence = NULL,
                            patent_incidence = NULL,
                            clinical_incidence = c(0, 1825),
                            severe_incidence = NULL,
                            prevalence = c(730, 3650),
                            hypnozoite_prevalence = NULL
){

  parent_formals <- names(formals())
  parent_formals <- parent_formals[which(parent_formals != "parameters")]
  outputs <- parent_formals[!unlist(lapply(parent_formals, function(x){is.null(get(x))}))]

  for (output in outputs) {
    parameters[[paste0(output, "_rendering_min_ages")]] <- get(output)[-length(get(output))]
    parameters[[paste0(output, "_rendering_max_ages")]] <- get(output)[-1]-1
  }

  parameters
}

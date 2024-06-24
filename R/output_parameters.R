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
#' @param ica age breaks for average acquired clinical immunity; default = NULL
#' @param icm age breaks for average maternal clinical immunity; default = NULL
#' @param iva age breaks for average acquired severe immunity; default = NULL
#' @param ivm age breaks for average maternal severe immunity; default = NULL
#' @param id age breaks for average immunity to detectability; default = NULL
#' @param ib age breaks for average blood immunity; default = NULL
#' @param hypnozoite_batches age breaks for average hypnozoite prevalence outputs (p.v only); default = NULL
#' @param n_with_hypnozoites age breaks for number of individuals with hypnozoites (p.v only); default = NULL
#' @export
#'
set_epi_outputs <- function(parameters,
                            age_group = NULL,
                            incidence = NULL,
                            patent_incidence = NULL,
                            clinical_incidence = NULL,
                            severe_incidence = NULL,
                            prevalence = NULL,
                            ica = NULL,
                            icm = NULL,
                            iva = NULL,
                            ivm = NULL,
                            id = NULL,
                            idm = NULL,
                            ib = NULL,
                            hypnozoites = NULL,
                            n_with_hypnozoites = NULL
){
  
  input <- list(
    age_group = age_group,
    incidence = incidence,
    patent_incidence = patent_incidence,
    clinical_incidence = clinical_incidence,
    severe_incidence = severe_incidence,
    prevalence = prevalence,
    ica = ica,
    icm = icm,
    iva = iva,
    ivm = ivm,
    id = id,
    idm = idm,
    ib = ib,
    hypnozoites = hypnozoites,
    n_with_hypnozoites = n_with_hypnozoites
  )
  input <- input[!sapply(input, is.null)]
  
  for(i in seq_along(input)){
    name <- names(input)[i]
    ages <- input[[i]]
    min_ages <- ages[-length(ages)]
    max_ages <- ages[-1] - 1
    parameters[[paste0(name, "_rendering_min_ages")]] <- min_ages
    parameters[[paste0(name, "_rendering_max_ages")]] <- max_ages
  }
  
  return(parameters)
}
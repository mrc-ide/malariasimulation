#' @title Parameterise age grouped output rendering
#'
#' @details this function produces discrete and contiguous age groups, inclusive of the lower
#' age limit and exclusive of the upper age limit: e.g., list(c(0, 10, 100), c(200, 250) will produce
#' three age groups: 0-9, 10-99 and 200-249 in days.
#' @param parameters the model parameters
#' @param age_group age breaks for population size outputs; default = NULL
#' @param incidence age breaks for incidence outputs (D+Tr+A); default = NULL
#' @param clinical_incidence age breaks for clinical incidence outputs (symptomatic); default = c(0, 1825)
#' @param severe_incidence age breaks for severe incidence outputs; default = NULL
#' @param prevalence age breaks for clinical prevalence outputs (pcr and lm detectable infections); default = c(730, 3650)
#' @param ica age breaks for average acquired clinical immunity; default = NULL
#' @param icm age breaks for average maternal clinical immunity; default = NULL
#' @param iva age breaks for average acquired severe immunity; default = NULL
#' @param ivm age breaks for average maternal severe immunity; default = NULL
#' @param id age breaks for average immunity to detectability; default = NULL
#' @param ib age breaks for average blood immunity; default = NULL
#' @export
#'
set_epi_outputs <- function(parameters,
                            age_group = NULL,
                            incidence = NULL,
                            clinical_incidence = NULL,
                            severe_incidence = NULL,
                            prevalence = NULL,
                            ica = NULL,
                            icm = NULL,
                            iva = NULL,
                            ivm = NULL,
                            id = NULL,
                            ib = NULL
){
  
  parent_formals <- names(formals())
  parent_formals <- parent_formals[which(parent_formals != "parameters")]
  outputs <- parent_formals[!unlist(lapply(parent_formals, function(x){is.null(get(x))}))]
  
  for (output in outputs) {
    if(!is.list(get(output))){stop("Each age input must be a list of vectors")}
    parameters[[paste0(output, "_rendering_min_ages")]] <- unlist(lapply(get(output), function(x){x[-length(x)]}))
    parameters[[paste0(output, "_rendering_max_ages")]] <- unlist(lapply(get(output), function(x){x[-1]-1}))
  }
  
  parameters
}

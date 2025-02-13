#' @title Parameterise age grouped output rendering
#'
#' @details this function produces contiguous age groups, inclusive of the lower
#' age limit and exclusive of the upper age limit: e.g., c(0, 10, 100) will produce
#' two age groups: 0-9 and 10-99 in days.
#' @param parameters the model parameters
#' @param age_group age breaks for population size outputs; default = NULL
#' @param incidence age breaks for incidence outputs (D+Tr+A); default = NULL
#' @param clinical_incidence age breaks for clinical incidence outputs (symptomatic); default = c(0, 1825)
#' @param severe_incidence age breaks for severe incidence outputs (pf only); default = NULL
#' @param prevalence age breaks for clinical prevalence outputs (pcr and lm detectable infections); default = c(730, 3650)
#' @param n_with_hypnozoites age breaks for individuals with hypnozoites outputs (pv only); default = NULL
#' @param ica age breaks for average acquired clinical immunity; default = NULL
#' @param icm age breaks for average maternal clinical immunity; default = NULL
#' @param iva age breaks for average acquired severe immunity (pf only); default = NULL
#' @param ivm age breaks for average maternal severe immunity (pf only); default = NULL
#' @param id age breaks for average immunity to detectability (pf only); default = NULL
#' @param ib age breaks for average blood immunity (pf only); default = NULL
#' @param iaa age breaks for average acquired antiparasite immunity (pv only); default = NULL
#' @param iam age breaks for average maternal antiparasite immunity (pv only); default = NULL
#' @param hypnozoites age breaks for average hypnozoite batches (pv only); default = NULL
#' @export
#'
set_epi_outputs <- function(parameters,
                            age_group = NULL,
                            incidence = NULL,
                            clinical_incidence = NULL,
                            severe_incidence = NULL,
                            prevalence = NULL,
                            n_with_hypnozoites = NULL,
                            ica = NULL,
                            icm = NULL,
                            iva = NULL,
                            ivm = NULL,
                            id = NULL,
                            ib = NULL,
                            iaa = NULL,
                            iam = NULL,
                            hypnozoites = NULL
){
  
  input <- list(
    age_group = age_group,
    incidence = incidence,
    clinical_incidence = clinical_incidence,
    severe_incidence = severe_incidence,
    prevalence = prevalence,
    n_with_hypnozoites = n_with_hypnozoites,
    ica = ica,
    icm = icm,
    iva = iva,
    ivm = ivm,
    id = id,
    ib = ib,
    iaa = iaa,
    iam = iam,
    hypnozoites = hypnozoites
  )
  
  if(parameters$parasite == "falciparum"){
    if(!any(lapply(list(n_with_hypnozoites, iaa, iam, hypnozoites), is.null))){
      stop("Some selected epi outputs do not apply to this parasite.")
    }
    
  } else if (parameters$parasite == "vivax"){
    if(!any(lapply(list(severe_incidence, iva, ivm, id, ib), is.null))){
      stop("Some selected epi outputs do not apply to this parasite.")
    }
  }
  
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

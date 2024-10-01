#' Parasite parameters
#'
#' Parasite-specific parameters for P. falciparum and vivax for use in malariasimulation
#'
#' @format ## `parasite_parameters`
#' A 
#'
#' @source <https://www.nature.com/articles/ncomms4136>, <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6097992/>
"parasite_parameters"

#' Parameter draws (P. falciparum)
#'
#' A list of 1000 draws from the joint posterior fit from
#'
#' @format ## `parameter_draws_pf`
#' A list of lists of length 1000, each level contains a list of drawn parameters
#'
#' @source <https://www.nature.com/articles/ncomms4136>
"parameter_draws_pf"

#' Parameter draws (P. vivax)
#'
#' 1000 draws from the joint posterior fit from
#'
#' @format ## `parameter_draws_pv`
#' A list of lists of length 1000, each level contains a list of drawn parameters
#'
#' @source <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6097992/>
"parameter_draws_pv"
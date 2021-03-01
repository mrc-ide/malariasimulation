#' @title Parameterise an RTS,S strategy
#' @param parameters a list of parameters to modify
#' @param timesteps a vector of timesteps for each round of vaccinations
#' @param coverages the coverage for each round of vaccinations
#' @param min_ages for the target population, inclusive (in timesteps)
#' @param max_ages for the target population, inclusive (in timesteps)
#' @param boosters the timesteps (following the initial vaccination) at which booster vaccinations are administered
#' @param booster_coverage the proportion of the vaccinated population who will
#' receive each booster vaccine
#' @export
set_rtss <- function(
  parameters,
  timesteps,
  coverages,
  min_ages,
  max_ages,
  boosters,
  booster_coverage
  ) {
  if (length(min_ages) != length(max_ages)) {
    stop('min and max ages do not align')
  }
  if (length(booster_coverage) != length(boosters)) {
    stop('booster and booster_coverage does not align')
  }
  parameters$rtss <- TRUE
  parameters$rtss_timesteps <- timesteps
  parameters$rtss_coverages <- coverages
  parameters$rtss_min_ages <- min_ages
  parameters$rtss_max_ages <- max_ages
  parameters$rtss_boosters <- boosters
  parameters$rtss_booster_coverage <- booster_coverage
  parameters
}

#' @title Parameterise an TBV strategy
#' @param parameters a list of parameters to modify
#' @param timesteps a vector of timesteps for each round of vaccinations
#' @param coverages the coverage for each round of vaccinations
#' @param ages for each round (in years)
#' vaccine
#' @export
set_tbv <- function(
  parameters,
  timesteps,
  coverages,
  ages
  ) {
  parameters$tbv <- TRUE
  parameters$tbv_timesteps <- timesteps
  parameters$tbv_coverages <- coverages
  parameters$tbv_ages <- ages
  parameters
}

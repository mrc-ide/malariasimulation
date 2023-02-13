#' @title Parameterise an RTS,S epi strategy
#'
#' @description distribute RTS,S vaccine when an individual becomes a certain
#' age. Efficacy will take effect after the last dose
#'
#' @param parameters a list of parameters to modify
#' @param coverages a vector of coverages for the starter doses
#' @param timesteps a vector of timesteps associated with coverages
#' @param age for the target population, (in timesteps)
#' @param min_wait the minimum acceptable time since the last vaccination (in
#' timesteps); When seasonal_boosters = TRUE, this represents the minimum time
#' between an individual receiving the final dose and the first booster. When using
#' both set_mass_rtss and set_rtss_epi, this represents the minimum time between
#' an individual being vaccinated under one scheme and vaccinated under another.
#' @param boosters the timesteps (following the final dose) at which booster vaccinations are administered
#' @param booster_coverage the proportion of the vaccinated population who will
#' receive each booster vaccine
#' @param seasonal_boosters logical, if TRUE the first booster timestep is
#' relative to the start of the year, otherwise they are relative to the last dose
#' @export
set_rtss_epi <- function(
  parameters,
  coverages,
  timesteps,
  age,
  min_wait,
  boosters,
  booster_coverage,
  seasonal_boosters = FALSE
  ) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  if(length(coverages) != length(timesteps)){
    stop("coverages and timesteps must align")
  }
  stopifnot(min_wait >= 0)
  stopifnot(age >= 0)
  stopifnot(is.logical(seasonal_boosters))
  if (seasonal_boosters) {
    if(boosters[[1]] < 0) {
      boosters <- boosters + 365
    }
  }
  stopifnot((length(boosters) == 0) || all(boosters > 0))
  stopifnot((length(booster_coverage)) == 0 || all(booster_coverage >= 0 & booster_coverage <= 1))
  if (length(booster_coverage) != length(boosters)) {
    stop('booster and booster_coverage does not align')
  }
  parameters$rtss <- TRUE
  parameters$rtss_epi_coverages <- coverages
  parameters$rtss_epi_timesteps <- timesteps
  parameters$rtss_epi_age <- age
  parameters$rtss_epi_boosters <- boosters
  parameters$rtss_epi_min_wait <- min_wait
  parameters$rtss_epi_booster_coverage <- booster_coverage
  parameters$rtss_epi_seasonal_boosters <- seasonal_boosters
  parameters
}

#' @title Parameterise an RTS,S mass distribution strategy
#'
#' @description distribute RTS,S vaccine to a population in an age range.
#' Efficacy will take effect after the last dose
#'
#' @param parameters a list of parameters to modify
#' @param timesteps a vector of timesteps for each round of vaccinations
#' @param coverages the coverage for each round of vaccinations
#' @param min_wait the minimum acceptable time since the last vaccination (in timesteps);
#' When using both set_mass_rtss and set_rtss_epi, this represents the minimum
#' time between an individual being vaccinated under one scheme and vaccinated under another.
#' @param min_ages for the target population, inclusive (in timesteps)
#' @param max_ages for the target population, inclusive (in timesteps)
#' @param boosters the timesteps (following the initial vaccination) at which booster vaccinations are administered
#' @param booster_coverage the proportion of the vaccinated population who will
#' receive each booster vaccine
#' @export
set_mass_rtss <- function(
  parameters,
  timesteps,
  coverages,
  min_ages,
  max_ages,
  min_wait,
  boosters,
  booster_coverage
  ) {
  stopifnot(all(timesteps > 1))
  stopifnot(min_wait >= 0)
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  stopifnot(all(min_ages >= 0 & max_ages >= 0))
  stopifnot(all(boosters > 0))
  stopifnot(all(booster_coverage >= 0 & booster_coverage <= 1))
  if (length(min_ages) != length(max_ages)) {
    stop('min and max ages do not align')
  }
  if (length(booster_coverage) != length(boosters)) {
    stop('booster and booster_coverage does not align')
  }
  parameters$rtss <- TRUE
  parameters$rtss_mass_timesteps <- timesteps
  parameters$rtss_mass_coverages <- coverages
  parameters$rtss_mass_min_ages <- min_ages
  parameters$rtss_mass_max_ages <- max_ages
  parameters$rtss_mass_min_wait <- min_wait
  parameters$rtss_mass_boosters <- boosters
  parameters$rtss_mass_booster_coverage <- booster_coverage
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
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  parameters$tbv <- TRUE
  parameters$tbv_timesteps <- timesteps
  parameters$tbv_coverages <- coverages
  parameters$tbv_ages <- ages
  parameters
}

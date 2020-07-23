#' @title Parameterise an RTS,S strategy
#' @param parameters a list of parameters to modify
#' @param start the timestep to start
#' @param end the last timestep for the intervention
#' @param frequency the number of timsteps between rounds
#' @param min_ages for the target population, inclusive (in timesteps)
#' @param max_ages for the target population, inclusive (in timesteps)
#' @param boosters the timesteps (following the initial vaccination) at which booster vaccinations are administered
#' @param coverage the proportion of the target population who will receive the
#' vaccine
#' @param booster_coverage the proportion of the vaccinated population who will
#' receive each booster vaccine
#' @export
set_rtss <- function(
  parameters,
  start,
  end,
  frequency,
  min_ages,
  max_ages,
  boosters,
  coverage,
  booster_coverage
  ) {
  if (length(min_ages) != length(max_ages)) {
    stop('min and max ages do not align')
  }
  if (start >= end) {
    stop('end must be strictly greater than start')
  }
  if (length(booster_coverage) != length(boosters)) {
    stop('booster and booster_coverage does not align')
  }
  parameters$rtss <- TRUE
  parameters$rtss_start <- start
  parameters$rtss_end <- end
  parameters$rtss_frequency <- frequency
  parameters$rtss_min_ages <- min_ages
  parameters$rtss_max_ages <- max_ages
  parameters$rtss_boosters <- boosters
  parameters$rtss_coverage <- coverage
  parameters$rtss_booster_coverage <- booster_coverage
  parameters
}

#' @title Parameterise an TBV strategy
#' @param parameters a list of parameters to modify
#' @param start the timestep to start
#' @param end the last timestep for the intervention
#' @param frequency the number of timsteps between rounds
#' @param ages for each round (in years)
#' @param coverage the proportion of the target population who will recieve the
#' vaccine
#' @export
set_tbv <- function(
  parameters,
  start,
  end,
  frequency,
  ages,
  coverage
  ) {
  parameters$tbv <- TRUE
  parameters$tbv_start <- start
  parameters$tbv_end <- end
  parameters$tbv_frequency <- frequency
  parameters$tbv_ages <- ages
  parameters$tbv_coverage <- coverage
  parameters
}

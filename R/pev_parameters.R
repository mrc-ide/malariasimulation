#' @title create a PEV profile
#' @description creates a data structure for holding pre-erythrocytic vaccine
#' profile parameters. Parameters are validated on creation.
#' @param vmax maximum efficacy of the vaccine
#' @param alpha shape parameter for the vaccine efficacy model
#' @param beta scale parameter for the vaccine efficacy model
#' @param cs peak parameters for the antibody model (vector of mean and std. dev)
#' @param rho delay parameters for the antibody model (vector of mean and std. dev)
#' @param ds delay parameters for the antibody model, short-term waning (vector of mean and std. dev)
#' @param dl delay parameters for the antibody model, long-term waning (vector of mean and std. dev)
#' @export
create_pev_profile <- function(vmax, alpha, beta, cs, rho, ds, dl) {
  allargs <- c(vmax, alpha, beta, cs, rho, ds, dl)
  stopifnot(all(is.numeric(allargs)))
  stopifnot(all(allargs > 0))
  stopifnot(length(cs) == 2)
  stopifnot(length(rho) == 2)
  stopifnot(length(ds) == 2)
  stopifnot(length(dl) == 2)
  list(
    vmax = vmax,
    alpha = alpha,
    beta = beta,
    cs = cs,
    rho = rho,
    ds = ds,
    dl = dl
  )
}

#' @title RTS,S vaccine profile
#' @description Parameters for a primary dose of RTS,S for use with the
#' set_mass_pev and set_pev_epi functions
#' @export
rtss_profile <- create_pev_profile(
  vmax = 0.93,
  alpha = 0.74,
  beta = 99.4,
  cs = c(6.37008, 0.35),
  rho = c(2.37832, 1.00813),
  ds = c(3.74502, 0.341185), # (White MT et al. 2015 Lancet ID)
  dl = c(6.30365, 0.396515) # (White MT et al. 2015 Lancet ID)
)

#' @title RTS,S booster vaccine profile
#' @description Parameters for a booster dose of RTS,S for use with the
#' set_mass_pev and set_pev_epi functions
#' @export
rtss_booster_profile <- create_pev_profile(
  vmax = 0.93,
  alpha = 0.74,
  beta = 99.4,
  cs = c(5.56277, 0.35),
  rho = c(1.03431, 1.02735),
  ds = c(3.74502, 0.341185), # (White MT et al. 2015 Lancet ID)
  dl = c(6.30365, 0.396515) # (White MT et al. 2015 Lancet ID)
)

#' @title Parameterise a pre-erythrocytic vaccine with an EPI strategy
#'
#' @description distribute vaccine when an individual becomes a certain
#' age. Efficacy will take effect after the last dose
#'
#' @param parameters a list of parameters to modify
#' @param profile primary vaccine profile of type PEVProfile
#' @param coverages a vector of coverages for the primary doses
#' @param timesteps a vector of timesteps associated with coverages
#' @param age the age when an individual will receive the first dose of the
#' vaccine (in timesteps)
#' @param min_wait the minimum acceptable time since the last vaccination (in
#' timesteps); When seasonal_boosters = TRUE, this represents the minimum time
#' between an individual receiving the final dose and the first booster. When using
#' both set_mass_pev and set_pev_epi, this represents the minimum time between
#' an individual being vaccinated under one scheme and vaccinated under another.
#' @param booster_timestep the timesteps (following the final dose) at which booster vaccinations are administered
#' @param booster_coverage the proportion of the vaccinated population relative to the last vaccination (whether a previous booster or the primary series)
#' @param booster_timed_coverage a time varying proportion of the vaccinated population relative to the last vaccination (whether a previous booster or the primary series), set in time with `booster_coverage_timestep`
#' @param booster_timed_coverage_timestep a vector of timesteps to change the time varying coverage specified in `booster_timed_coverage`
#' @param booster_profile list of booster vaccine profiles, of type
#' PEVProfile, for each timestep in booster_timeteps
#' @param seasonal_boosters logical, if TRUE the first booster timestep is
#' relative to the start of the year, otherwise they are relative to the last dose
#' @export
set_pev_epi <- function(
  parameters,
  profile,
  coverages,
  timesteps,
  age,
  min_wait,
  booster_timestep,
  booster_coverage,
  booster_timed_coverage = NULL,
  booster_timed_coverage_timestep = NULL,
  booster_profile,
  seasonal_boosters = FALSE
  ) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))

  # Check that the primary timing parameters make sense
  if(length(coverages) != length(timesteps)){
    stop("coverages and timesteps must align")
  }

  # Check that seasonal booster parameters make sense
  stopifnot(min_wait >= 0)
  stopifnot(age >= 0)
  stopifnot(is.logical(seasonal_boosters))
  if (seasonal_boosters) {
    if(booster_timestep[[1]] < 0) {
      booster_timestep <- booster_timestep + 365
    }
  }

  # Check that the booster timing parameters make sense
  stopifnot((length(booster_timestep) == 0) || all(booster_timestep > 0))
  stopifnot((length(booster_coverage)) == 0 || all(booster_coverage >= 0 & booster_coverage <= 1))
  if (!all(c(length(booster_coverage), length(booster_timestep), length(booster_profile)) == length(booster_timestep))) {
    stop('booster_timestep and booster_coverage and booster_profile does not align')
  }
  if (length(booster_timed_coverage) != length(booster_timed_coverage_timestep)) {
    stop("booster_coverage_timestep must be the same length as booster_coverage")
  }

  # Index the new vaccine profiles
  profile_list <- c(list(profile), booster_profile)
  profile_indices <- seq_along(profile_list) + length(parameters$pev_profiles)
  parameters$pev_profiles <- c(parameters$pev_profiles, profile_list)

  # Set the EPI strategy
  parameters$pev <- TRUE
  parameters$pev_epi_coverages <- coverages
  parameters$pev_epi_timesteps <- timesteps
  parameters$pev_epi_age <- age
  parameters$pev_epi_booster_timestep <- booster_timestep
  parameters$pev_epi_min_wait <- min_wait
  parameters$pev_epi_booster_coverage <- booster_coverage
  parameters$pev_epi_timed_booster_coverage <- booster_timed_coverage
  parameters$pev_epi_timed_booster_coverage_timestep <- booster_timed_coverage_timestep
  parameters$pev_epi_booster_coverage <- booster_coverage
  parameters$pev_epi_profile_indices <- profile_indices
  parameters$pev_epi_seasonal_boosters <- seasonal_boosters
  parameters
}

#' @title Parameterise a vaccine mass distribution strategy
#'
#' @description distribute pre-erythrocytic vaccine to a population in an age range.
#' Efficacy will take effect after the last dose
#'
#' @param parameters a list of parameters to modify
#' @param profile primary vaccine profile of type PEVProfile
#' @param timesteps a vector of timesteps for each round of vaccinations
#' @param coverages the coverage for each round of vaccinations
#' @param min_wait the minimum acceptable time since the last vaccination (in timesteps);
#' When using both set_mass_pev and set_pev_epi, this represents the minimum
#' time between an individual being vaccinated under one scheme and vaccinated under another.
#' @param min_ages for the target population, inclusive (in timesteps)
#' @param max_ages for the target population, inclusive (in timesteps)
#' @param booster_timestep the timesteps (following the initial vaccination) at which booster vaccinations are administered
#' @param booster_coverage the proportion of the vaccinated population relative to the last vaccination (whether a previous booster or the primary series)
#' @param booster_timed_coverage a time varying proportion of the vaccinated population relative to the last vaccination (whether a previous booster or the primary series), set in time with `booster_coverage_timestep`
#' @param booster_timed_coverage_timestep a vector of timesteps to change the time varying coverage specified in `booster_timed_coverage`
#' @param booster_profile list of booster vaccine profiles, of type
#' PEVProfile, for each timestep in booster_timeteps
#' @export
set_mass_pev <- function(
  parameters,
  profile,
  timesteps,
  coverages,
  min_ages,
  max_ages,
  min_wait,
  booster_timestep,
  booster_coverage,
  booster_timed_coverage = NULL,
  booster_timed_coverage_timestep = NULL,
  booster_profile
  ) {
  stopifnot(all(timesteps >= 1))
  stopifnot(min_wait >= 0)
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  stopifnot(all(min_ages >= 0 & max_ages >= 0))
  stopifnot(all(booster_timestep > 0))
  stopifnot(all(booster_coverage >= 0 & booster_coverage <= 1))
  if (length(min_ages) != length(max_ages)) {
    stop('min and max ages do not align')
  }
  if (!all(c(length(booster_coverage), length(booster_timestep), length(booster_profile)) == length(booster_timestep))) {
    stop('booster_timestep, booster_coverage and booster_profile does not align')
  }
  if (length(booster_timed_coverage) != length(booster_timed_coverage_timestep)) {
    stop("booster_coverage_timestep must be the same length as booster_coverage")
  }

  # Index the new vaccine profiles
  profile_list <- c(list(profile), booster_profile)
  profile_indices <- seq_along(profile_list) + length(parameters$pev_profiles)
  parameters$pev_profiles <- c(parameters$pev_profiles, profile_list)

  # Set the mass vaccination strategy
  parameters$pev <- TRUE
  parameters$mass_pev_timesteps <- timesteps
  parameters$mass_pev_coverages <- coverages
  parameters$mass_pev_min_ages <- min_ages
  parameters$mass_pev_max_ages <- max_ages
  parameters$mass_pev_min_wait <- min_wait
  parameters$mass_pev_booster_timestep <- booster_timestep
  parameters$mass_pev_booster_coverage <- booster_coverage
  parameters$mass_pev_booster_timed_coverage <- booster_timed_coverage
  parameters$mass_pev_booster_timed_coverage_timestep <- booster_timed_coverage_timestep
  parameters$mass_pev_profile_indices <- profile_indices
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

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
#' set_mass_pev and set_pev_epi functions (White MT et al. 2015 Lancet ID)
#' @export
rtss_profile <- create_pev_profile(
  vmax = 0.93,
  alpha = 0.74,
  beta = 99.4,
  cs = c(6.37008, 0.35),
  rho = c(2.37832, 1.00813),
  ds = c(3.74502, 0.341185),
  dl = c(6.30365, 0.396515)
)

#' @title R21 vaccine profile
#' @description Parameters for a primary dose of R21 for use with the
#' set_mass_pev and set_pev_epi functions (Schmit + Topazian et al. 2022 Lancet ID)
#' @export
r21_profile <- create_pev_profile(
  vmax = 0.87,
  alpha = 0.91,
  beta = 471,
  cs = c(9.3199677, 0.8387902),
  rho = c(0.8071676, 0.6010363),
  ds = c(3.7996007, 0.1618982),
  dl = c(6.2820200, 0.4549185)
)


#' @title RTS,S booster vaccine profile
#' @description Parameters for a booster dose of RTS,S for use with the
#' set_mass_pev and set_pev_epi functions (White MT et al. 2015 Lancet ID)
#' @export
rtss_booster_profile <- create_pev_profile(
  vmax = 0.93,
  alpha = 0.74,
  beta = 99.4,
  cs = c(5.56277, 0.35),
  rho = c(1.03431, 1.02735),
  ds = c(3.74502, 0.341185),
  dl = c(6.30365, 0.396515)
)

#' @title R21 booster vaccine profile
#' @description Parameters for a booster dose of R21 for use with the
#' set_mass_pev and set_pev_epi functions (Schmit + Topazian et al. 2022 Lancet ID)
#' @export
r21_booster_profile <- create_pev_profile(
  vmax = 0.87,
  alpha = 0.91,
  beta = 471,
  cs = c(9.2372858, 0.7188541),
  rho = c(0.07140337, 0.54175154),
  ds = c(3.7996007, 0.1618982), 
  dl = c(6.2820200, 0.4549185)
)

#' @title Parameterise a pre-erythrocytic vaccine with an EPI strategy
#'
#' @description distribute vaccine when an individual becomes a certain
#' age. Efficacy will take effect after the last dose
#'
#' @param parameters a list of parameters to modify
#' @param profile a list of details for the vaccine profile, create with `create_pev_profile`
#' @param coverages a vector of coverages for the primary doses
#' @param timesteps a vector of timesteps for each change in coverage
#' @param age the age when an individual will receive the first dose of the
#' vaccine (in timesteps)
#' @param min_wait the minimum acceptable time since the last vaccination (in
#' timesteps); When seasonal_boosters = TRUE, this represents the minimum time
#' between an individual receiving the final dose and the first booster. When using
#' both set_mass_pev and set_pev_epi, this represents the minimum time between
#' an individual being vaccinated under one scheme and vaccinated under another.
#' @param booster_spacing the timesteps (following the final primary dose) at which booster vaccinations are administered
#' @param booster_coverage a matrix of coverages (timesteps x boosters) specifying the proportion the previously vaccinated population to continue receiving booster doses. The rows of the matrix must be the same size as `timesteps`. The columns of the matrix must be the same size as `booster_spacing`.
#' @param booster_profile list of lists representing each booster profile, the outer list must be the same length as `booster_spacing`. Create vaccine profiles with `create_pev_profile`
#' @param seasonal_boosters logical, if TRUE the first booster timestep is
#' relative to the start of the year, otherwise they are relative to the last primary dose
#' @export
set_pev_epi <- function(
  parameters,
  profile,
  coverages,
  timesteps,
  age,
  min_wait,
  booster_spacing,
  booster_coverage,
  booster_profile,
  seasonal_boosters = FALSE
  ) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  stopifnot(is.matrix(booster_coverage))

  # Check that the primary timing parameters make sense
  if(length(coverages) != length(timesteps)){
    stop("coverages and timesteps must align")
  }

  # Check that booster_spacing are monotonically increasing
  if (length(booster_spacing) > 1) {
    if (!all(diff(booster_spacing) > 0)) {
      stop('booster_spacing must be monotonically increasing')
    }
  }

  # Check that seasonal booster parameters make sense
  stopifnot(min_wait >= 0)
  stopifnot(age >= 0)
  stopifnot(is.logical(seasonal_boosters))
  if (seasonal_boosters) {
    if(booster_spacing[[1]] < 0) {
      booster_spacing <- booster_spacing + 365
    }
  }

  # Check that the booster timing parameters make sense
  stopifnot((length(booster_spacing) == 0) || all(booster_spacing > 0))
  stopifnot((length(booster_coverage)) == 0 || all(booster_coverage >= 0 & booster_coverage <= 1))
  if (!all(c(ncol(booster_coverage), length(booster_profile)) == length(booster_spacing))) {
    stop('booster_spacing, booster_coverage and booster_profile do not align')
  }
  # Check that booster_coverage and timesteps align
  if (length(booster_coverage) > 0) {
    if (nrow(booster_coverage) != length(timesteps)) {
      stop('booster_coverage and timesteps do not align')
    }
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
  parameters$pev_epi_booster_spacing <- booster_spacing
  parameters$pev_epi_min_wait <- min_wait
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
#' @param profile a list of details for the vaccine profile, create with `create_pev_profile`
#' @param timesteps a vector of timesteps for each round of vaccinations
#' @param coverages the coverage for each round of vaccinations
#' @param min_wait the minimum acceptable time since the last vaccination (in timesteps);
#' When using both set_mass_pev and set_pev_epi, this represents the minimum
#' time between an individual being vaccinated under one scheme and vaccinated under another.
#' @param min_ages for the target population, inclusive (in timesteps)
#' @param max_ages for the target population, inclusive (in timesteps)
#' @param booster_spacing the timesteps (following the final primary dose) at which booster vaccinations are administered
#' @param booster_coverage a matrix of coverages (timesteps x boosters) specifying the proportion the previously vaccinated population to continue receiving booster doses. The rows of the matrix must be the same size as `timesteps`. The columns of the matrix must be the same size as `booster_spacing`.
#' @param booster_profile list of lists representing each booster profile, the outer list must be the same length as `booster_spacing`. Create vaccine profiles with `create_pev_profile`
#' @export
set_mass_pev <- function(
  parameters,
  profile,
  timesteps,
  coverages,
  min_ages,
  max_ages,
  min_wait,
  booster_spacing,
  booster_coverage,
  booster_profile
  ) {
  stopifnot(all(timesteps >= 1))
  stopifnot(min_wait >= 0)
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  stopifnot(all(min_ages >= 0 & max_ages >= 0))
  stopifnot(all(booster_spacing > 0))
  stopifnot(all(booster_coverage >= 0 & booster_coverage <= 1))
  if (length(min_ages) != length(max_ages)) {
    stop('min and max ages do not align')
  }

  # Check that booster_spacing are monotonically increasing
  if (length(booster_spacing) > 1) {
    if (!all(diff(booster_spacing) > 0)) {
      stop('booster_spacing must be monotonically increasing')
    }
  }

  stopifnot((length(booster_coverage)) == 0 || all(booster_coverage >= 0 & booster_coverage <= 1))
  if (!all(c(ncol(booster_coverage), length(booster_profile)) == length(booster_spacing))) {
    stop('booster_spacing, booster_coverage and booster_profile do not align')
  }
  # Check that booster_coverage and timesteps align
  if (length(booster_coverage) > 0) {
    if (nrow(booster_coverage) != length(timesteps)) {
      stop('booster_coverage and timesteps do not align')
    }
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
  parameters$mass_pev_booster_spacing <- booster_spacing
  parameters$mass_pev_booster_coverage <- booster_coverage
  parameters$mass_pev_profile_indices <- profile_indices
  parameters
}

#' @title Vaccine profile
#' @description A data structure for holding vaccine profile parameters.
#' Parameters are validated on creation.
#' @importFrom R6 R6Class
#' @export
VaccineProfile <- R6::R6Class(
  'VaccineProfile',
  public = list(
    #' @field vmax maximum efficacy of the vaccine
    vmax,
    #' @field alpha shape parameter for the vaccine efficacy model
    alpha,
    #' @field beta scale parameter for the vaccine efficacy model
    beta,
    #' @field cs peak parameters for the antibody model (vector of mean and std. dev)
    cs,
    #' @field rho delay parameters for the antibody model (vector of mean and std. dev)
    rho,
    #' @field ds delay parameters for the antibody model, short-term waning (vector of mean and std. dev)
    ds,
    #' @field dl delay parameters for the antibody model, long-term waning (vector of mean and std. dev)
    dl,

    #' @description create a vaccine profile
    #' @param vmax immutable value for vmax
    #' @param alpha immutable value for alpha
    #' @param beta immutable value for beta
    #' @param cs immutable value for cs
    #' @param rho immutable value for rho
    #' @param ds immutable value for ds
    #' @param dl immutable value for dl
    initialize = function(vmax, alpha, beta, cs, rho, ds, dl) {
      allargs <- c(vmax, alpha, beta, cs, rho, ds, dl)
      stopifnot(all(is.numeric(allargs)))
      stopifnot(all(allargs > 0))
      stopifnot(length(cs) == 2)
      stopifnot(length(rho) == 2)
      stopifnot(length(ds) == 2)
      stopifnot(length(dl) == 2)
      self$vmax = vmax
      self$alpha = alpha
      self$beta = beta
      self$cs = cs
      self$rho = rho
      self$ds = ds
      self$dl = dl
    }
  )
)

#' @export
rtss_profile <- VaccineProfile$new(
  vmax = 0.93
  alpha = 0.74
  beta = 99.4
  cs = c(6.37008, 0.35)
  rho = c(2.37832, 1.00813)
  ds = c(3.74502, 0.341185) # (White MT et al. 2015 Lancet ID)
  dl = c(6.30365, 0.396515) # (White MT et al. 2015 Lancet ID)
)

#' @export
rtss_booster_profile <- VaccineProfile$new(
  vmax = 0.93
  alpha = 0.74
  beta = 99.4
  cs = c(5.56277, 0.35)
  rho = c(1.03431, 1.02735)
  ds = c(3.74502, 0.341185) # (White MT et al. 2015 Lancet ID)
  dl = c(6.30365, 0.396515) # (White MT et al. 2015 Lancet ID)
)

#' @title Parameterise a vaccine epi strategy
#'
#' @description distribute vaccine when an individual becomes a certain
#' age. Efficacy will take effect after the last dose
#'
#' @param parameters a list of parameters to modify
#' @param profile primary vaccine profile of type VaccineProfile
#' @param coverages a vector of coverages for the primary doses
#' @param timesteps a vector of timesteps associated with coverages
#' @param age for the target population, (in timesteps)
#' @param min_wait the minimum acceptable time since the last vaccination (in
#' timesteps); When seasonal_boosters = TRUE, this represents the minimum time
#' between an individual receiving the final dose and the first booster. When using
#' both set_mass_vaccine and set_vaccine_epi, this represents the minimum time between
#' an individual being vaccinated under one scheme and vaccinated under another.
#' @param booster_timestep the timesteps (following the final dose) at which booster vaccinations are administered
#' @param booster_coverage the proportion of the vaccinated population who will
#' receive each booster vaccine
#' @param booster_profile list of booster vaccine profiles, of type
#' VaccineProfile, for each timestep in booster_timeteps
#' @param seasonal_boosters logical, if TRUE the first booster timestep is
#' relative to the start of the year, otherwise they are relative to the last dose
#' @export
set_vaccine_epi <- function(
  parameters,
  profile,
  coverages,
  timesteps,
  age,
  min_wait,
  booster_timestep,
  booster_coverage,
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
    if(boosters[[1]] < 0) {
      boosters <- boosters + 365
    }
  }

  # Check that the booster timing parameters make sense
  stopifnot((length(booster_timestep) == 0) || all(boosters > 0))
  stopifnot((length(booster_coverage)) == 0 || all(booster_coverage >= 0 & booster_coverage <= 1))
  if (!all(c(length(booster_coverage), length(booster_timestep), length(booster_profile)) == length(booster_timestep))) {
    stop('booster_timestep and booster_coverage and booster_profile does not align')
  }

  # Check that vaccine profiles are well formed
  stopifnot(inherits(profile, 'VaccineProfile'))
  for (p in booster_profile) {
    stopifnot(inherits(p, 'VaccineProfile'))
  }

  parameters$vaccines <- TRUE
  parameters$vaccine_epi_profile <- profile
  parameters$vaccine_epi_coverages <- coverages
  parameters$vaccine_epi_timesteps <- timesteps
  parameters$vaccine_epi_age <- age
  parameters$vaccine_epi_booster_timestep <- booster_timestep
  parameters$vaccine_epi_min_wait <- min_wait
  parameters$vaccine_epi_booster_coverage <- booster_coverage
  parameters$vaccine_epi_booster_profile <- booster_profile
  parameters$vaccine_epi_seasonal_boosters <- seasonal_boosters
  parameters
}

#' @title Parameterise a vaccine mass distribution strategy
#'
#' @description distribute vaccine to a population in an age range.
#' Efficacy will take effect after the last dose
#'
#' @param parameters a list of parameters to modify
#' @param profile primary vaccine profile of type VaccineProfile
#' @param timesteps a vector of timesteps for each round of vaccinations
#' @param coverages the coverage for each round of vaccinations
#' @param min_wait the minimum acceptable time since the last vaccination (in timesteps);
#' When using both set_mass_rtss and set_rtss_epi, this represents the minimum
#' time between an individual being vaccinated under one scheme and vaccinated under another.
#' @param min_ages for the target population, inclusive (in timesteps)
#' @param max_ages for the target population, inclusive (in timesteps)
#' @param booster_timestep the timesteps (following the initial vaccination) at which booster vaccinations are administered
#' @param booster_coverage the proportion of the vaccinated population who will
#' receive each booster vaccine
#' @param booster_profile list of booster vaccine profiles, of type
#' VaccineProfile, for each timestep in booster_timeteps
#' @export
set_mass_vaccination <- function(
  parameters,
  profile,
  timesteps,
  coverages,
  min_ages,
  max_ages,
  min_wait,
  booster_timestep,
  booster_coverage
  booster_profile,
  ) {
  stopifnot(all(timesteps > 1))
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

  # Check that vaccine profiles are well formed
  stopifnot(inherits(profile, 'VaccineProfile'))
  for (p in booster_profile) {
    stopifnot(inherits(p, 'VaccineProfile'))
  }
  parameters$vaccines <- TRUE
  parameters$mass_vaccine_timesteps <- timesteps
  parameters$mass_vaccine_coverages <- coverages
  parameters$mass_vaccine_min_ages <- min_ages
  parameters$mass_vaccine_max_ages <- max_ages
  parameters$mass_vaccine_min_wait <- min_wait
  parameters$mass_vaccine_booster_timestep <- booster_timestep
  parameters$mass_vaccine_booster_coverage <- booster_coverage
  parameters$mass_vaccine_booster_profile <- booster_profile
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

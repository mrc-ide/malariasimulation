#' @title Calculate equilibrium solution for vector counts
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
#' @param foim equilibrium foim
#' @param m (optional) the total number of female adult mosquitos
initial_mosquito_counts <- function(parameters, foim = 0, m = NULL) {
  if (is.null(m)) {
    m <- parameters$total_M
  }
  omega <- calculate_omega(parameters)
  n_E <- 2 * omega * parameters$mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_L <- 2 * parameters$mum * parameters$dl * (
    1. + parameters$dpl * parameters$mup
  ) * m

  n_P <- 2 * parameters$dpl * parameters$mum * m

  n_Sm <- m * parameters$mum / (foim + parameters$mum)

  incubation_survival <- exp(-parameters$mum * parameters$dem)

  n_Pm <- m * foim / (foim + parameters$mum) * (
    1. - incubation_survival
  )

  n_Im <- m * foim / (foim + parameters$mum) * incubation_survival

  c(n_E, n_L, n_P, n_Sm, n_Pm, n_Im)
}

#' @title Calculate omega value
#' @description useful value for calculating equilibrium solutions for vectors
#' taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
calculate_omega <- function(parameters) {
  sub_omega <- parameters$gamma * parameters$ml / parameters$me - (
    parameters$del / parameters$dl
  ) + (
    (parameters$gamma - 1) * parameters$ml * parameters$del
  )

  -.5 * sub_omega + sqrt(
    .25 * sub_omega**2 +
      .5 * parameters$gamma * parameters$beta * parameters$ml * parameters$del /
      (parameters$me * parameters$mum * parameters$dl * (
        1. + parameters$dpl * parameters$mup
      ))
  )
}

#' @title Calculate the vector carrying capacity
#' @description taken from
#' "Modelling the impact of vector control interventions on Anopheles gambiae
#' population dynamics"
#' @param parameters model parameters
calculate_carrying_capacity <- function(parameters) {
  m <- parameters$total_M
  omega <- calculate_omega(parameters)

  m * 2 * parameters$dl * parameters$mum * (
    1. + parameters$dpl * parameters$mup
  ) * parameters$gamma * (omega + 1) / (
    omega / (parameters$ml * parameters$del) - (
      1. / (parameters$ml * parameters$dl)
    ) - 1.
  )
}

#' @title Calculate the mean rainfall throughout the year
#' @param parameters model parameters
calculate_R_bar <- function(parameters) {
  mean(vnapply(1:365, function(t) rainfall(
		t,
		parameters$days_per_timestep,
    parameters$g0,
    parameters$g,
    parameters$h
	)))
}

#' @title Calculate equilibrium total_M from parameters
#'
#' @param parameters to work from
#' @param EIR equilibrium to use, bites per person per year
equilibrium_total_M <- function(parameters, EIR) {
  if (EIR == 0) {
    return(0)
  }
  if (parameters$init_foim == 0) {
    stop('init_foim must be > 0 to calculate a non-zero equilibrium total_M')
  }
  total_daily_eir <- EIR * parameters$human_population / 365
  lifetime <- parameters$init_foim * exp(-parameters$mum * parameters$dem) / (
    parameters$init_foim + parameters$mum
  )
  total_daily_eir / sum(
    parameters$variety_proportions * parameters$blood_meal_rates * parameters$Q0 * lifetime
  )
}

#' @title Calculate the yearly offset (in timesteps) for the peak mosquito
#' season
#'
#' @param parameters to work from
#' @export
peak_season_offset <- function(parameters) {
  K0 <- calculate_carrying_capacity(parameters)
  R_bar <- calculate_R_bar(parameters)
  which.max(vnapply(seq(365), function(t) {
    carrying_capacity(
      t,
      parameters$model_seasonality,
      parameters$days_per_timestep,
      parameters$g0,
      parameters$g,
      parameters$h,
      K0,
      R_bar
    )
  }))[[1]]
}

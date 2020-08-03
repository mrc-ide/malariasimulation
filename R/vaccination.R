calculate_rtss_antibodies <- function(
  t,
  cs,
  rho,
  ds,
  dl,
  parameters
  ) {
  cs * (
    rho * exp(-t * log(2) / ds) + (
      1 - rho
    ) * exp(-t * log(2) / dl)
  )
}

calculate_rtss_efficacy <- function(antibodies, parameters) {
  parameters$rtss_vmax * (
    1 - (1 / (
      1 + (antibodies / parameters$rtss_beta) ** parameters$rtss_alpha
    ))
  )
}

calculate_tbv_antibodies <- function(t, parameters) {
  parameters$tbv_tau * (
    parameters$tbv_rho * exp(-t * log(2) / parameters$tbv_ds) + (
      1 - parameters$tbv_rho
    ) * exp(-t * log(2) / parameters$tbv_dl)
  )
}

calculate_TRA <- function(antibodies, parameters) {
  numerator <- (
    antibodies / parameters$tbv_tra_mu
  ) ** parameters$tbv_gamma1 
  numerator / (numerator + parameters$tbv_gamma2)
}

calculate_TBA <- function(z, mx, k) {
  scale <- 1 / (1 - (k / (k + mx)) ** k)
  offset <- (k / (k + mx)) ** k
  tra_transformation <- (k / (k + mx * (1 - z))) ** k
  scale * (tra_transformation - offset)
}

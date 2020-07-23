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

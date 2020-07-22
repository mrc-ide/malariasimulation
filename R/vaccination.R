calculate_rtss_antibodies <- function(t, boosted, parameters) {
  antibodies <- rep(NA, length(t))
  antibodies[!boosted] <- calculate_rtss_antibodies_from_params(
    t,
    parameters$cs_peak,
    parameters$rho_peak,
    parameters
  )
  antibodies[boosted] <- calculate_rtss_antibodies_from_params(
    t,
    parameters$cs_boost,
    parameters$rho_boost,
    parameters
  )
  antibodies
}

calculate_rtss_antibodies_from_params <- function(
  t,
  cs,
  rho,
  parameters
  ) {
  cs * (
    rho * exp(-t * log(2) / parameters$rtss_ds) + (
      1 - rho
    ) * exp(-t * log(2) / parameters$rtss_dl)
  )
}

calculate_rtss_efficacy <- function(antibodies, parameters) {
  parameters$v_max * (
    1 - (1 / (
      1 + (antibodies / parameters$rtss_beta) ** parameters$rtss_alpha
    ))
  )
}

#' @title Calculate the effect of TBV on infectivity
#' @description Returns a vector of human infectivity towards mosquitoes
#' accounting for the reduction in transmission due to vaccination
#'
#' @param timestep current timestep
#' @param infectivity a vector of raw infectivities
#' @param variables the available variables
#' @param parameters model parameters
account_for_tbv <- function(
  timestep,
  infectivity,
  variables,
  parameters
) {
  time_vaccinated <- variables$tbv_vaccinated$get_values()
  vaccinated <- which(time_vaccinated != -1)
  affected_states <- c('U', 'A', 'D', 'Tr')
  mx <- parameters[c('tbv_mu', 'tbv_ma', 'tbv_md', 'tbv_mt')]
  for (i in seq_along(affected_states)) {
    in_state <- variables$state$get_index_of(affected_states[[i]])$to_vector()
    vaccinated_in_state <- intersect(vaccinated, in_state)
    antibodies <- calculate_tbv_antibodies(
      timestep - time_vaccinated[vaccinated_in_state],
      parameters$tbv_tau,
      parameters$tbv_rho,
      parameters$tbv_ds,
      parameters$tbv_dl
    )
    tra <- calculate_TRA(
      parameters$tbv_tra_mu,
      parameters$tbv_gamma1,
      parameters$tbv_gamma2,
      antibodies
    )
    tba <- calculate_TBA(
      mx[[i]],
      parameters$tbv_k,
      tra
    )
    infectivity[vaccinated_in_state] <- infectivity[vaccinated_in_state] * (
      1 - tba
    )
  }
  infectivity
}

calculate_tbv_antibodies <- function(t, tau, rho, ds, dl) {
  tau * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl))
}

calculate_TRA <- function(mu, gamma1, gamma2, antibodies) {
  numerator <- (antibodies / mu)^gamma1
  numerator / (numerator + gamma2)
}

calculate_TBA <- function(mx, k, tra) {
  offset <- (k / (k + mx)) ^ k;
  scale <- 1. / (1. - offset);
  tra_transformation <- (k / (k + mx * (1 - tra))) ^ k;
  scale * (tra_transformation - offset)
}


#' @title return the probability of being bitten given vector controls
#' @param individuals a list of available individuals
#' @param variables a list of available variables
#' @param Im the infectious mosquito state
#' @param api the simulation api
#' @param parameters model parameters
prob_bitten <- function(individuals, variables, species, api, parameters) {
  if (!(parameters$bednets || parameters$spraying)) {
    return(rep(1, api$get_parameters()$human_population))
  }

  timestep <- api$get_timestep()

  if (parameters$bednets) {
    phi_bednets <- parameters$phi_bednets[[species]]
    net_time <- api$get_variable(individuals$human, variables$net_time)
    sn <- prob_survives_bednets(timestep - net_time, species, parameters)
    sn[net_time == -1] <- 1
  } else {
    phi_bednets <- 0
    sn <- 1
  }

  if (parameters$spraying) {
    phi_spraying <- parameters$phi_spraying[[species]]
    spray_time <- api$get_variable(individuals$human, variables$spray_time)
    rs <- prob_spraying_repels(
      timestep - spray_time,
      parameters$rs[[species]],
      parameters$gammas
    )
    rs_ <- 1 - rs # the complement of rs
    ss <- prob_survives_spraying(
      spray_time,
      parameters$endophily[[species]],
      parameters$gammas
    )
    ss[spray_time == -1] <- 1
  } else {
    phi_spraying <- 0
    rs_ <- 1
    ss <- 1
  }
  
  return(
    1 - phi_spraying +
      phi_bednets * rs_ * sn +
      (phi_spraying - phi_bednets) * rs_
  )
}

prob_spraying_repels <- function(t, rs0, gammas) {
  rs0 * vector_control_decay(t, gammas)
}

prob_survives_spraying <- function(t, endophily, gammas) {
  ds <- endophily * vector_control_decay(t, gammas)
  1 - ds
}

prob_survives_bednets <- function(t, species, parameters) {
  rnm <- parameters$rnm[[species]]
  gamman <- parameters$gamman
  rn <- (parameters$rn[[species]] - rnm) * vector_control_decay(t, gamman) + rnm
  dn <- parameters$dn0[[species]] * vector_control_decay(t, gamman)
  1 - rn - dn
}

vector_control_decay <- function(t, gamma) {
  exp(-t * gamma)
}

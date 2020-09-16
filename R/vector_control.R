
#' @title return the probability of being bitten given vector controls
#' @param individuals a list of available individuals
#' @param variables a list of available variables
#' @param species the species to calculate for
#' @param api the simulation api
#' @param parameters model parameters
prob_bitten <- function(individuals, variables, species, api, parameters) {
  if (!(parameters$bednets || parameters$spraying)) {
    n <- parameters$human_population
    return(
      list(
        prob_bitten_survives = rep(1, n),
        prob_bitten = rep(1, n),
        prob_repelled = rep(0, n)
      )
    )
  }

  timestep <- api$get_timestep()

  if (parameters$bednets) {
    phi_bednets <- parameters$phi_bednets[[species]]
    net_time <- api$get_variable(individuals$human, variables$net_time)
    since_net <- timestep - net_time
    rn <- prob_repelled_bednets(since_net, species, parameters)
    sn <- prob_survives_bednets(rn, since_net, species, parameters)
    sn[net_time == -1] <- 1
    rn[net_time == -1] <- 0
  } else {
    phi_bednets <- 0
    sn <- 1
    rn <- 0
  }

  if (parameters$spraying) {
    phi_spraying <- parameters$phi_spraying[[species]]
    spray_time <- api$get_variable(individuals$human, variables$spray_time)
    since_spray <- timestep - spray_time
    rs <- prob_spraying_repels(
      since_spray,
      parameters$rs[[species]],
      parameters$gammas
    )
    rs_ <- 1 - rs # the complement of rs
    ss <- prob_survives_spraying(
      since_spray,
      parameters$endophily[[species]],
      parameters$gammas
    )
    ss[spray_time == -1] <- 1
  } else {
    phi_spraying <- 0
    rs <- 0
    rs_ <- 1
    ss <- 1
  }
  
  list(
    prob_bitten_survives = (
      1 - phi_spraying +
      phi_bednets * rs_ * sn * ss +
      (phi_spraying - phi_bednets) * rs_ * ss
    ),
    prob_bitten = (
      1 - phi_spraying +
      phi_bednets * rs_ * sn +
      (phi_spraying - phi_bednets) * rs_
    ),
    prob_repelled = (
      phi_bednets * rs_ * rn +
      phi_spraying * rs
    )
  )
}

prob_spraying_repels <- function(t, rs0, gammas) {
  rs0 * vector_control_decay(t, gammas)
}

prob_survives_spraying <- function(t, endophily, gammas) {
  ds <- endophily * vector_control_decay(t, gammas)
  1 - ds
}

prob_repelled_bednets <- function(t, species, parameters) {
  rnm <- parameters$rnm[[species]]
  gamman <- parameters$gamman
  (parameters$rn[[species]] - rnm) * vector_control_decay(t, gamman) + rnm
}

prob_survives_bednets <- function(rn, t, species, parameters) {
  dn <- parameters$dn0[[species]] * vector_control_decay(t, parameters$gamman)
  1 - rn - dn
}

vector_control_decay <- function(t, gamma) {
  exp(-t * gamma)
}

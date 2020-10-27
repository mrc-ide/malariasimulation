
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
    unused <- net_time == -1
    api$render('net_usage', sum(!unused) / length(unused))
    sn[unused] <- 1
    rn[unused] <- 0
  } else {
    phi_bednets <- 0
    sn <- 1
    rn <- 0
  }

  if (parameters$spraying) {
    phi_indoors <- parameters$phi_indoors[[species]]
    spray_time <- api$get_variable(individuals$human, variables$spray_time)
    since_spray <- timestep - spray_time
    unused <- spray_time == -1
    rs <- prob_spraying_repels(
      since_spray,
      parameters$rs[[species]],
      parameters$gammas
    )
    rs[unused] <- 0
    rs_comp <- 1 - rs
    ss <- prob_survives_spraying(
      since_spray,
      parameters$endophily[[species]],
      parameters$gammas
    )
    ss[unused] <- 1
  } else {
    phi_indoors <- 0
    rs <- 0
    rs_comp <- 1
    ss <- 1
  }

  
  list(
    prob_bitten_survives = (
      1 - phi_indoors +
      phi_bednets * rs_comp * sn * ss +
      (phi_indoors - phi_bednets) * rs_comp * ss
    ),
    prob_bitten = (
      1 - phi_indoors +
      phi_bednets * rs_comp * sn +
      (phi_indoors - phi_bednets) * rs_comp
    ),
    prob_repelled = (
      phi_bednets * rs_comp * rn +
      phi_indoors * rs
    )
  )
}

indoor_spraying <- function(human, spray_time, parameters) {
  function(api) {
    timestep <- api$get_timestep()
    matches <- timestep == parameters$spraying_timesteps
    if (any(matches)) {
      target <- bernoulli(
        parameters$human_population,
        parameters$spraying_coverages[matches]
      )
      api$queue_variable_update(human, spray_time, timestep, target)
    }
  }
}

distribute_nets <- function(human, variables, parameters) {
  function(api) {
    timestep <- api$get_timestep()
    matches <- timestep == parameters$bednet_timesteps
    if (any(matches)) {
      target <- bernoulli(
        parameters$human_population,
        parameters$bednet_coverages[matches]
      )
      api$queue_variable_update(human, variables$net_time, timestep, target)
      api$queue_variable_update(
        human,
        variables$net_end_time,
        timestep + ceiling(
          -parameters$bednet_retention * log(runif(length(target)))
        ),
        target
      )
    }
  }
}

throw_away_nets <- function(human, variables) {
  function(api) {
    end_time <- api$get_variable(human, variables$net_end_time)
    target <- which(end_time == api$get_timestep())
    api$queue_variable_update(human, variables$net_time, -1, target)
    api$queue_variable_update(human, variables$net_end_time, -1, target)
  }
}

# =================
# Utility functions
# =================
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
  exp(-t / gamma)
}

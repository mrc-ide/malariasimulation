#' @title The probability of being bitten given vector controls
#' @param timestep current timestep
#' @param variables a list of available variables
#' @param species the species to calculate for
#' @param parameters model parameters
#' @noRd
prob_bitten <- function(
  timestep,
  variables,
  species,
  parameters
  ) {
  n <- parameters$human_population
  if (!(parameters$bednets || parameters$spraying || parameters$housing)) {
    return(
      list(
        prob_bitten_survives = rep(1, n),
        prob_bitten = rep(1, n),
        prob_repelled = rep(0, n)
      )
    )
  }

  if (parameters$bednets) {
    phi_bednets <- parameters$phi_bednets[[species]]
    net_time <- variables$net_time$get_values()
    since_net <- timestep - net_time
    matches <- match(net_time, parameters$bednet_timesteps)
    rn <- prob_repelled_bednets(matches, since_net, species, parameters)
    sn <- prob_survives_bednets(rn, matches, since_net, species, parameters)
    unused <- net_time == -1
    sn[unused] <- 1
    rn[unused] <- 0
  } else {
    phi_bednets <- 0
    sn <- 1
    rn <- 0
  }

  if (parameters$spraying) {
    phi_indoors <- parameters$phi_indoors[[species]]
    protected <- variables$spray_time$get_index_of(set=-1)$not(TRUE)
    spray_time <- variables$spray_time$get_values(protected)
    matches <- match(spray_time, parameters$spraying_timesteps)
    ls_theta <- parameters$spraying_ls_theta[matches, species]
    ls_gamma <- parameters$spraying_ls_gamma[matches, species]
    ks_theta <- parameters$spraying_ks_theta[matches, species]
    ks_gamma <- parameters$spraying_ks_gamma[matches, species]
    ms_theta <- parameters$spraying_ms_theta[matches, species]
    ms_gamma <- parameters$spraying_ms_gamma[matches, species]
    since_spray <- timestep - spray_time
    ls <- spraying_decay(since_spray, ls_theta, ls_gamma)
    ks <- parameters$k0 * spraying_decay(since_spray, ks_theta, ks_gamma)
    ms <- spraying_decay(since_spray, ms_theta, ms_gamma)
    js <- 1 - ls - ks
    ms_comp <- (1 - ms)
    ls_prime <- ls * ms_comp
    ks_prime <- ks * ms_comp
    js_prime <- js * ms_comp + ms
    protected_index <- protected$to_vector()
    rs <- rep(0, n)
    rs[protected_index] <- prob_spraying_repels(
      ls_prime,
      ks_prime,
      js_prime,
      parameters$k0
    )
    rs_comp <- 1 - rs
    ss <- rep(1, n)
    ss[protected_index] <- prob_survives_spraying(
      ks_prime,
      parameters$k0
    )
  } else {
    # phi_indoors <- 0 ## if housing on too this has to be phi_indoors
    rs <- 0
    rs_comp <- 1
    ss <- 1
  }
  
  if (parameters$housing) {
    house_time <- variables$house_time$get_values()
    since_house <- timestep - house_time
    matches <- match(house_time, parameters$house_timesteps)
    phi_housing <- parameters$phi_housing[matches, species]
    # rn_house <- parameters$rn_house[matches, species]
    # phi_housing <- parameters$phi_housing[[species]]
    rn_house <- prob_repelled_house(matches, since_house, species, parameters)
    sn <- 1 - rn
    unused <- house_time == -1
  } else {
    phi_housing <- 1
    rn_house <- 0
  }
  
  if (!(parameters$housing & parameters$spraying)) {
    phi_indoors <- 0 ## we want phi_indoors to be applied if housing is on
  }

  list(
    prob_bitten_survives = (
      1 - phi_indoors * phi_housing +
        (1 - rn_house) * (phi_bednets * phi_housing * rs_comp * sn * ss) +
        (1 - rn_house) * ((phi_indoors * phi_housing - phi_bednets * phi_housing) * rs_comp * ss)
    ),
    prob_bitten = (
      1 - phi_indoors * phi_housing +
        (1 - rn_house) * (phi_bednets * phi_housing * rs_comp * sn) +
        (1 - rn_house) * ((phi_indoors * phi_housing - phi_bednets * phi_housing) * rs_comp)
    ),
    prob_repelled = (
      rn_house * phi_indoors * phi_housing + 
        (1 - rn_house) * phi_bednets * phi_housing * rs_comp * rn + 
        (1 - rn_house) * phi_indoors * phi_housing * rs
    )
  )
}

#' @title Simulate housing improvements
#' @description simulates improved housing so that harder for vectors to get indoors
#' from `set_housing` and correlation parameters from
#' `get_correlation_parameters`
#'
#' @param variables list of variables in the model
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
housing_improvement <- function(variables, parameters, correlations) {
  function(timestep) {
    matches <- timestep == parameters$house_timesteps
    if (any(matches)) {
      target <- which(sample_intervention(
        seq(parameters$human_population),
        'housing',
        parameters$house_coverages[matches],
        correlations
      ))
      house_time$queue_update(timestep, target)
    }
  }
}

#' @title Indoor spraying
#' @description models indoor residual spraying according to the strategy
#' from `set_spraying` and correlation parameters from
#' `get_correlation_parameters`
#'
#' @param spray_time the variable for the time of spraying
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
indoor_spraying <- function(spray_time, parameters, correlations) {
      function(timestep) {
    matches <- timestep == parameters$spraying_timesteps
    if (any(matches)) {
      target <- which(sample_intervention(
        seq(parameters$human_population),
        'spraying',
        parameters$spraying_coverages[matches],
        correlations
      ))
      spray_time$queue_update(timestep, target)
    }
  }
}

#' @title Distribute nets
#' @description distributes nets to individuals according to the strategy
#' from `set_bednets` and correlation parameters from
#' `get_correlation_parameters`
#'
#' @param variables list of variables in the model
#' @param throw_away_net an event to trigger when the net will be removed
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
distribute_nets <- function(variables, throw_away_net, parameters, correlations) {
  function(timestep) {
    matches <- timestep == parameters$bednet_timesteps
    if (any(matches)) {
      target <- which(sample_intervention(
        seq(parameters$human_population),
        'bednets',
        parameters$bednet_coverages[matches],
        correlations
      ))
      variables$net_time$queue_update(timestep, target)
      throw_away_net$clear_schedule(target)
      throw_away_net$schedule(
        target,
        log_uniform(length(target), parameters$bednet_retention)
      )
    }
  }
}

throw_away_nets <- function(variables) {
  function(timestep, target) {
    variables$net_time$queue_update(-1, target)
  }
}

# =================
# Utility functions
# =================
prob_spraying_repels <- function(ls_prime, ks_prime, js_prime, k0) {
  (1 - ks_prime / k0) * (js_prime / (ls_prime + js_prime))
}

prob_survives_spraying <- function(ks_prime, k0) {
  ks_prime / k0
}

prob_repelled_bednets <- function(matches, dt, species, parameters) {
  rnm <- parameters$bednet_rnm[matches, species]
  gamman <- parameters$bednet_gamman[matches]
  (parameters$bednet_rn[matches, species] - rnm) * bednet_decay(dt, gamman) + rnm
}

prob_survives_bednets <- function(rn, matches, dt, species, parameters) {
  dn0 <- parameters$bednet_dn0[matches, species]
  dn <- dn0 * bednet_decay(dt, parameters$bednet_gamman[matches])
  1 - rn - dn
}

bednet_decay <- function(t, gamma) {
  exp(-t / gamma)
}

spraying_decay <- function(t, theta, gamma) {
  1 / (1 + exp(-(theta + gamma * t)))
}

prob_repelled_house <- function(t, matches, species, parameters) {
  parameters$rn_house[matches, species] * t/t    ## make this continual through time or we could have a decay as nets but very long lasting?
}

net_usage_renderer <- function(net_time, renderer) {
  function(t) {
    renderer$render(
      'n_use_net',
      net_time$get_index_of(-1)$not(TRUE)$size(),
      t
    )
  }
}

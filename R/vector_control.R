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
  if (!(parameters$bednets || parameters$spraying || parameters$spatial_emanator || parameters$spatial_emanator_outdoor)) {
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
    # adjust_se_check = 0
  } else {
    phi_bednets <- 0
    sn <- 1
    rn <- 0
    # adjust_se_check = 1
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
    spray_on = 1
    rs_comp <- 1 - rs
    ss <- rep(1, n)
    ss[protected_index] <- prob_survives_spraying(
      ks_prime,
      parameters$k0
    )
  } else {
    rs <- 0
    rs_comp <- 1
    ss <- 1
  }

  if (parameters$spatial_emanator) {
    
    phi_bednets <- parameters$phi_bednets[[species]]
    phi_indoors <- parameters$phi_indoors[[species]]
    
    ## Identify individuals who are protected by SEs
    protected <- variables$spatial_emanator_time$get_index_of(set=-1)$not(TRUE)
    
    ## Get time since SEs deployed
    spatial_emanator_time <- variables$spatial_emanator_time$get_values(protected)
    
    ## Get parameter values for deployed SEs
    matches <- match(spatial_emanator_time, parameters$spatial_emanator_timesteps)
    
    ## parameters indoor repellence (R), indoor killing before feeding (K), 
    ## indoor killing after feeding (M) and outdoor repellence (M)
    R0 <- parameters$spatial_emanator_R0[matches, species]
    K0 <- parameters$spatial_emanator_K0[matches, species]
    M0 <- parameters$spatial_emanator_M0[matches, species]
    delta0 <- parameters$spatial_emanator_delta0[matches, species]
    nR <- parameters$spatial_emanator_nR[matches, species]
    nK <- parameters$spatial_emanator_nK[matches, species]
    nM <- parameters$spatial_emanator_nM[matches, species]
    nd <- parameters$spatial_emanator_nd[matches, species]
    eps_R <- parameters$spatial_emanator_eps_R[matches, species]
    eps_K <- parameters$spatial_emanator_eps_K[matches, species]
    eps_M <- parameters$spatial_emanator_eps_M[matches, species]
    eps_d <- parameters$spatial_emanator_eps_d[matches, species]
    T_ref <- parameters$spatial_emanator_T_ref[matches, species]
    
    since_spatial_emanator <- timestep - spatial_emanator_time
    
    ## Translate input parameters into time dependent trends
    mu_R <- hill_eps(R0, T_ref, nR, since_spatial_emanator, eps_R)
    mu_K <- hill_eps(K0, T_ref, nK, since_spatial_emanator, eps_K)
    mu_M <- hill_eps(M0, T_ref, nM, since_spatial_emanator, eps_M)
    mu_d <- hill_eps(delta0, T_ref, nd, since_spatial_emanator, eps_d)
    
    ## Translate time-dependent trends into model probabilities
    
    # Set whole population values
    protected_index <- protected$to_vector()
    rse <- kse <- bse <- mse <- sse <- rep(0, n)
    
    ## At the first decision point a mosquito may be repelled (r), may be killed (k) or may bite (b)
    rse[protected_index] <- mu_d + (1-mu_d) * mu_R
    kse[protected_index] <- (1-mu_d) * mu_K
    bse <- 1 - rse - kse
    
    ## Given that a mosquito has bitten, the mosquito may still be killed or survive
    mse[protected_index] <- mu_M
    sse <- 1 - mse

    spatial_emanator_on <- 1
    
    } else {
      spatial_emanator_on <- 0
      rse <- 0 ## repellence
      kse <- 0 ## killed before feeding
      bse <- 1 ## bites
      mse <- 0 ## killed after biting
      sse <- 1 ## survives after biting
  }
  
  if ((!parameters$spatial_emanator & !parameters$spraying)) {
    phi_indoors <- 0 ## we want phi_indoors to be applied if spatial emanators is on
  }

  # if (parameters$spatial_emanator_outdoor) {
  #   phi_bednets <- parameters$phi_bednets[[species]]
  #   phi_indoors <- parameters$phi_indoors[[species]]
  # 
  #   protected <- variables$spatial_emanator_outdoor_time$get_index_of(set=-1)$not(TRUE)
  #   spatial_emanator_outdoor_time <- variables$spatial_emanator_outdoor_time$get_values(protected)
  #   matches <- match(spatial_emanator_outdoor_time, parameters$spatial_emanator_outdoor_timesteps)
  # 
  #   rse_out_theta <- parameters$spatial_emanator_out_theta[matches, species]
  #   rse_out_gamma <- parameters$spatial_emanator_out_gamma[matches, species]
  # 
  #   dse_out_theta <- parameters$spatial_emanator_mort_out_theta[matches, species]
  #   dse_out_gamma <- parameters$spatial_emanator_mort_out_gamma[matches, species]
  # 
  #   since_spatial_emanator_outdoor <- timestep - spatial_emanator_outdoor_time
  # 
  #   rse_out_1 <- spraying_decay(since_spatial_emanator_outdoor, rse_out_theta, rse_out_gamma)
  #   dse_out_1 <- spraying_decay(since_spatial_emanator_outdoor, dse_out_theta, dse_out_gamma)
  # 
  #   rse_out_temp <- 1 - rse_out_1
  # 
  #   protected_index <- protected$to_vector()
  #   rse_out <- rep(0, n)
  #   rse_out[protected_index] <- rse_out_1
  #   dse_out <- rep(0, n)
  #   dse_out[protected_index] <- rse_out_temp * dse_out_1
  # 
  #   dfse_out <- dse_out * 0.04 ## fraction of mosquitoes having fed then died
  #   rse_out_comp <- rep(0, n)
  # 
  #   sse_out <- 1 - rse_out - dse_out
  #   rse_out_comp <- 1 - rse_out
  # 
  # } else {
  #   rse_out <- 0
  #   dfse_out <- 0
  #   dse_out <- 0
  #   rse_out_comp <- 1
  #   sse_out <- 1
  # }
  
  # list(
  #   prob_bitten_survives = (
  #     (1 - phi_indoors) * sse_out +
  #     phi_bednets * rs_comp * sn * ss * sse_in +
  #     (phi_indoors - phi_bednets) * rs_comp * ss * sse_in
  #   ),
  #   prob_bitten = (
  #     (1 - phi_indoors) * (sse_out + dfse_out) +
  #     phi_bednets * rs_comp * sn * (sse_in + dfse_in) +
  #     (phi_indoors - phi_bednets) * rs_comp * (sse_in + dfse_in)
  #   ),
  #   prob_repelled = (
  #     phi_bednets * rs_comp * rn * rse_in_comp +
  #       phi_indoors * repel_inside +         ## will go to 0 if IRS off
  #       (1 - phi_indoors) * rse_out                    ## 0 if outdoor emanators off
  #   )
  # )

  list(
    prob_bitten_survives = (
      
      # Outdoor biting: assumes outdoor mosquitos bite
      (1 - phi_indoors) +
        
        # Indoor biting: assumes mosquitos are not repelled by spraying (rs_comp) and
        # are not repelled or killed prior to biting by SE (bse)
        # and survive resting on sprayed walls following feeding (ss)
        # and survive SE killing after feeding (sse)
        (phi_indoors - phi_bednets) * rs_comp * bse * ss * sse +
        
        # Bed biting: assumes mosquitos are not repelled by spraying (rs_comp) and
        # are not repelled or killed before biting by SE (bse)
        # and survive nets (sn) and survive spraying (ss) and survives SE (sse)
        phi_bednets * rs_comp * bse * sn * ss * sse
      
    ),
    
    prob_bitten = (
      
      # Outdoor biting: assumes outdoor mosquitos bite
      (1 - phi_indoors) +
        
        # Indoor biting: assumes mosquitos are not repelled by spraying (rs_comp) and 
        # are not repelled or killed by SEs prior to biting (bse)
        (phi_indoors - phi_bednets) * rs_comp * bse +
        
        # Bet biting: assumes mosquitos are not repelled by spraying (rs_comp) and 
        # are not repelled or killed by SEs prior to biting (bse)
        # and are not repelled or killed by nets (sn)
        phi_bednets * rs_comp * bse * sn
      
      ## IRS parameters: 
      ## rs + rs_comp = 1, ds + ss = 1, rs + ds + ss != 1
      
      ## Similarly, se parameters: 
      ## rse + rse_comp = 1, rse + dse + x (survive first stage to bite) = 1, dfse + sse = 1, rse + ds + ss != 1
      
      ## If we have rse + dse + dfse + sse = 1, we need to group them carefully
      ## rse + rse_comp = 1, where rse_comp = dse + dfse + sse
      ## If a mosquito makes it through to biting, the probability is: 1 - rse - dse = dfse + sse = rse_comp - dse
      ## Any of these should be fine to use in the conditional probability that you survive, given that you have bitten.
      ## e.g. sse_in/(1 - rse - dse) = sse_in/(dfse + sse) = sse_in/(rse_comp - dse)
    ),
    
    prob_repelled = (
      
      # Outdoor biting: assume no repellence
      # (1 - phi_indoors) +
        
        # Indoor biting: may be affected by spraying (s) and SE (se)
        # We assume these are independent events
        phi_indoors * (1 - (1 - rs) * (1 - rse)) +
        
        # Bet biting: assumes that mosquitos have not been repelled by spraying (rs_comp)
        # or repelled of killed by SE (bse) prior to being repelled by nets (rn)
        phi_bednets * rs_comp * bse * rn
      
    )
  )
  
}


#' @title Spatial Emanators Outdoors
#' @description models outdoor use of spatial emanators according to the strategy
#' from `set_spatial_emanator_outdoor` and correlation parameters from
#' `get_correlation_parameters`
#'
#' @param spatial_emanator_outdoor_time the variable for the time of spatial emanator deployment
#' @param renderer model rendering object
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
spatial_emanator_outdoor <- function(spatial_emanator_outdoor_time, renderer, parameters, correlations) {
  renderer$set_default('n_spatial_emanator_outdoor', 0)
  function(timestep) {
    matches <- timestep == parameters$spatial_emanator_outdoor_timesteps
    if (any(matches)) {
      target <- which(sample_intervention(
        seq(parameters$human_population),
        'spatial_emanator_outdoor',
        parameters$spatial_emanator_outdoor_coverages[matches],
        correlations
      ))
      spatial_emanator_outdoor_time$queue_update(timestep, target)
      renderer$render('n_spatial_emanator_outdoor', length(target), timestep)
    }
  }
}

#' @title Spatial Emanators
#' @description models indoor use of spatial emanators according to the strategy
#' from `set_spatial_emanator` and correlation parameters from
#' `get_correlation_parameters`
#'
#' @param spatial_emanator_time the variable for the time of spatial emanator deployment
#' @param renderer model rendering object
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
spatial_emanator <- function(spatial_emanator_time, renderer, parameters, correlations) {
  renderer$set_default('n_spatial_emanator', 0)
  function(timestep) {
    matches <- timestep == parameters$spatial_emanator_timesteps
    if (any(matches)) {
      target <- which(sample_intervention(
        seq(parameters$human_population),
        'spatial_emanator',
        parameters$spatial_emanator_coverages[matches],
        correlations
      ))
      spatial_emanator_time$queue_update(timestep, target)
      renderer$render('n_spatial_emanator', length(target), timestep)
    }
  }
}

#' @title Indoor spraying
#' @description models indoor residual spraying according to the strategy
#' from `set_spraying` and correlation parameters from
#' `get_correlation_parameters`
#'
#' @param spray_time the variable for the time of spraying
#' @param renderer model rendering object
#' @param parameters the model parameters
#' @param correlations correlation parameters
#' @noRd
indoor_spraying <- function(spray_time, renderer, parameters, correlations) {
  renderer$set_default('n_spray', 0)
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
      renderer$render('n_spray', length(target), timestep)
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
# prob_spatial_emanator_decay <- function(matches, dt, species, rse_in, parameters){
#   dse_in_theta = parameters$dse_in_theta[matches, species]
#   dse_in_gamma = parameters$dse_in_gamma[matches, species]
#   (1 - rse_in) * spraying_decay(dt,dse_in_theta, dse_in_gamma)
# }


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

hill_eps <- function(f0, T_ref, n, t, eps){
  rho <- 1/eps - 1
  return(f0 / (1 + rho * (t / T_ref)^n))
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

#' @title Calculate the effect of TBV on infectivity
#' @description Returns a vector of human infectivity towards mosquitoes
#' accounting for the reduction in transmission due to vaccination
#'
#' @param timestep current timestep
#' @param infectivity a vector of raw infectivities
#' @param variables the available variables
#' @param parameters model parameters
#' @noRd
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
    #JDC:Do we need the equivalent of 'vaccine_efficacy' used in RTSS?
    #JDC: What is the equivalent of 'source_vector' in this case?
    vaccine_times <- variables$tbv$get_values(source_vector) #JDC. source_vector(or equivalent) missing at present
    vaccination_index <- source_vector[vaccine_times > -1] #JDC source_vector(or equivalent) missing at present
    antibodies <- calculate_tbv_antibodies(
      timestep - time_vaccinated[vaccinated_in_state],
      #parameters$tbv_iiv
      variables$tbv_iiv$get_values(vaccinated_index) #JDC: IIV parameters added
      #parameters$tbv_tau,
      #parameters$tbv_rho,
      #parameters$tbv_ds,
      #parameters$tbv_dl
    )
    tra <- calculate_TRA(
      #parameters$tbv_tra_mu,
      #parameters$tbv_gamma1,
      #parameters$tbv_gamma2,
      parameters$tbv_hill,
      parameters$tbv_EC50,
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


#' @title Distribute TBV vaccine
#'
#' @param timestep current timestep
#' @param variables the available variables
#' @param events the available events
#' @param parameters model parameters
#' @param correlations model correlations
#' @param renderer object for model outputs
#' @noRd
create_tbv_listener <- function(variables, events, parameters, correlations, renderer) {
  function(timestep) {
    time_index = which(parameters$tbv_timesteps == timestep)
    target <- which(trunc(get_age(
      variables$birth$get_values(),
      timestep
    ) / 365) %in% parameters$tbv_ages)
    to_vaccinate <- which(sample_intervention(
      target,
      'tbv',
      parameters$tbv_coverages[[time_index]],
      correlations
    ))
    renderer$render('n_vaccinated_tbv', length(to_vaccinate), timestep)
    if (length(to_vaccinate) > 0) {
      variables$tbv_vaccinated$queue_update(
        timestep,
        to_vaccinate
      )
    }
    if (time_index < length(parameters$tbv_timesteps)) {
      events$tbv_vaccination$schedule(
        parameters$tbv_timesteps[[time_index + 1]] - timestep
      )
    }
  }
}

calculate_tbv_antibodies <- function(t, tbv_iiv){ #Note: tbv_iiv has length 2
  #tau * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl))

  #params
  TVCL <- 0.0051
  TVV1 <- 2.57
  TVV2 <- 0.944
  TVQ <- 0.172
  TVV3 <- 1.47
  TVQ2 <- 0.0782
  WT <- 70 #body weight (fixed for now, will have to adjust)
  AMT <- WT*10

  WTCL <- (WT/70)**0.75
  WTV <- (WT/70)**1

  CL <- TVCL*exp(tbv_iiv[1])*WTCL #IIV
  V1 <- TVV1*exp(tbv_iiv[2])*WTV #IIV
  V2 <- TVV2*WTV
  Q <- TVQ*WTCL
  V3 <- TVV3*WTV
  Q2 <- TVQ2*WTCL

  #Rate consts.
  k10 <- 24*CL/V1
  k12 <- 24*Q/V1
  k21 <- 24*Q/V2
  k13 <- 24*Q2/V1
  k31 <- 24*Q2/V3

  amt <- 700 # dose will eventually be age-dependent

  mt <- matrix(c(-(k12+k10+k13), k21, k31, k12, -k21, 0, k13, 0, -k31), nrow = 3, byrow = T)
  p <- eigen(mt)$vectors #These are in the right order
  lambda <- eigen(mt)$values
  y0 <- c(amt,0,0) 
  zz0 <- solve(p)%*%y0
  z <- c(zz0[1]*exp(lambda[1] * max( t ,.1)), #time fudge! (avoid t=0)
         zz0[2]*exp(lambda[2] * max( t ,.1)),
         zz0[3]*exp(lambda[3] * max( t ,.1)))
  zz <-  p[1,1]*z[1] + p[1,2]*z[2] + p[1,3]*z[3]# p%*%z #
  zz/V1
}

calculate_TRA <- function(hill, EC50, antibodies) { 
  #numerator <- (antibodies / mu)^gamma1
  #numerator / (numerator + gamma2)
  (antibodies^hill)/((antibodies^hill) + (EC50^hill))
}

calculate_TBA <- function(mx, k, tra) { #same as before
  offset <- (k / (k + mx)) ^ k
  scale <- 1. / (1. - offset)
  tra_transformation <- (k / (k + mx * (1 - tra))) ^ k
  scale * (tra_transformation - offset)
}
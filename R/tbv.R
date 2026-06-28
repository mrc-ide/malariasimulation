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
  #time_vaccinated <- variables$tbv_vaccinated$get_values()
  #vaccinated <- which(time_vaccinated != -1)
  # much faster than which
  vaccinated <- variables$tbv_vaccinated$get_index_of(set=-1)$not(TRUE) #JDC: added a 'TRUE' here too, as was producing warning message
  affected_states <- c('U', 'A', 'D', 'Tr')
  #mx <- parameters[c('tbv_mu', 'tbv_ma', 'tbv_md', 'tbv_mt')]
  for (i in seq_along(affected_states)) {
    #in_state <- variables$state$get_index_of(affected_states[[i]])$to_vector()
    in_state <- variables$state$get_index_of(affected_states[[i]]) # keep as a bitset?
    #vaccinated_in_state <- intersect(vaccinated, in_state)
    vaccinated_in_state <- in_state$and(vaccinated) # much faster than intersect
    vaccine_times <- variables$tbv_vaccinated$get_values(vaccinated_in_state) # JDC: added back in (was in TB31F). I think this has been replaced with 'time_vaccinated'
    antibodies <- calculate_tbv_antibodies(
      timestep - vaccine_times, #time_vaccinated[vaccinated_in_state],
      variables$tbv_iiv1$get_values(vaccinated_in_state),
      variables$tbv_iiv2$get_values(vaccinated_in_state),
      variables$tbv_iiv3$get_values(vaccinated_in_state),
      variables$tbv_iiv4$get_values(vaccinated_in_state),
      parameters$tbv_adult_scaling,
      #parameters$tbv_tau,
      #parameters$tbv_rho,
      #parameters$tbv_ds,
      #parameters$tbv_dl,
      get_age(variables$birth$get_values(vaccinated_in_state), vaccine_times) #Will this select the relevant time?
      #variables$tbv_PK_zscore$get_values(vaccinated_in_state)
    )
    # tra <- calculate_TRA(
    #parameters$tbv_tra_mu,
    #parameters$tbv_gamma1,
    #parameters$tbv_gamma2,
    #  antibodies
    #)
    tba <- calculate_TBA(
      antibodies
      #parameters
      #mx[[i]],
      #parameters$tbv_k,
      #tra
    )
    vaccinated_vector <- vaccinated_in_state$to_vector() # JDC needs to be a vector!!
    infectivity[vaccinated_vector] <- infectivity[vaccinated_vector] * (
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
    to_vaccinate <- target[sample_intervention(
      target,
      'tbv',
      parameters$tbv_coverages[[time_index]],
      correlations
    )]
    renderer$render('n_vaccinated_tbv', length(to_vaccinate), timestep)
    if (length(to_vaccinate) > 0) {
      variables$tbv_vaccinated$queue_update(
        timestep,
        to_vaccinate
      )
    }
  }
}

#JDC: Add Age at vaccination
#calculate_tbv_antibodies <- function(t, tau, rho, ds, dl, adult_scaling){ # , age_at_vaccination. random variables, too!
calculate_tbv_antibodies <- function(t, tbv_iiv1, tbv_iiv2, tbv_iiv3, tbv_iiv4, adult_scaling, age_at_vaccination){ 
  #scaling <- 1
  #if(age_at_vaccination > 16*365){
  #  scaling = adult_scaling
  #}
  #peak3 <- 1558
  #peak4 <- 3648
  #tau <- 0.5*(peak3 + peak4) # or halfway???
  #use vnapply
  vnapply(
    #seq_along(tbv_iiv1),
    seq_along(age_at_vaccination),
    function(i) {
      #unpack params
      iiv1 <- tbv_iiv1[[i]] # rs
      iiv2 <- tbv_iiv2[[i]] # rl
      iiv3 <- tbv_iiv3[[i]] # rho
      iiv4 <- tbv_iiv4[[i]] # peak
      age_agent <- age_at_vaccination[[i]]
      tt <- t[[i]] #JDC: new addition- t now a vector? Perhaps because vaccination time may not be the same for everyone
      #adult scaling not a vector. Is this a problem
      rs <- exp(3.74 + 0.59*iiv1)
      rl <- exp(5.79 + 0.81*iiv2)
      rho1 <- invlogit(0.765 + 0.725*iiv3)
      rho2 <- invlogit(1.15 + 0.68*iiv3)
      peak1 <- exp(7.34 + 0.680*iiv4) # (3rd dose)
      peak2 <- exp(8.202 + 0.677*iiv4) # (4th dose)
      peak <- peak1
      rho <- rho1
      ttt <- tt
      if(tt > 365){
        peak <- peak2
        rho <- rho2
        ttt <- tt - 365
      }
      scaling <- 1
      if(age_agent > 16*365){
        scaling <- adult_scaling
      }
      peak * scaling * (rho * exp(-ttt/rs) + (1 - rho) * exp(-ttt/rl))
      # tau * adult_scaling * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl))
      #tau * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl))
    }
  )
}

#calculate_TRA <- function(mu, gamma1, gamma2, antibodies) {
#JDC: could remove
calculate_TRA <- function(antibodies) {
  #numerator <- (antibodies / mu)^gamma1
  #numerator / (numerator + gamma2)
  antibodies
}

calculate_TBA <- function(antibodies){ #}, parameters) {
  # offset <- (k / (k + mx)) ^ k;
  # scale <- 1. / (1. - offset);
  # tra_transformation <- (k / (k + mx * (1 - tra))) ^ k;
  # scale * (tra_transformation - offset)
  #JDC: Pfs230 model: 
  hill <- 6.2 #ish
  c50 <- 134 #ish
  Vmax <- 0.921 #ish
  Vmax*((antibodies^hill)/((c50^hill)+(antibodies^hill)))
}
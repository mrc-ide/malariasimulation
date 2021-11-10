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
    #vaccinated <- which(time_vaccinated != -1)
    # much faster than which
    vaccinated <- variables$tbv_vaccinated$get_index_of(set=-1)$not(TRUE) #JDC: added a 'TRUE' here too, as was producing warning message
    # NOTE: this requires changing tbv_vaccinated to an IntegerVariable 
    affected_states <- c('U', 'A', 'D', 'Tr')
    mx <- parameters[c('tbv_mu', 'tbv_ma', 'tbv_md', 'tbv_mt')]
    for (i in seq_along(affected_states)) {
      #in_state <- variables$state$get_index_of(affected_states[[i]])$to_vector()
      in_state <- variables$state$get_index_of(affected_states[[i]]) #keep as a Bitset
      #vaccinated_in_state <- intersect(vaccinated, in_state)
      vaccinated_in_state <- in_state$and(vaccinated) # much faster than intersect
      #age_vector <- get_age(variables$birth$get_values(vaccinated_in_state), timestep) ?
      #vaccine_times <- variables$tbv$get_values(source_vector)
      vaccine_times <- variables$tbv_vaccinated$get_values(vaccinated_in_state)
      antibodies <- calculate_tbv_antibodies(
        timestep - vaccine_times,
        variables$tbv_iiva$get_values(vaccinated_in_state),
        variables$tbv_iivb$get_values(vaccinated_in_state),
        variables$tbv_PK_sx$get_values(vaccinated_in_state),
        variables$tbv_PK_zscore$get_values(vaccinated_in_state), #then add age
        get_age(variables$birth$get_values(vaccinated_in_state), vaccine_times) #Will this select the relevant time?
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
    vaccinated_vector <- vaccinated_in_state$to_vector()
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

calculate_tbv_antibodies <- function(t, tbv_iiva, tbv_iivb, tbv_PK_sx, tbv_PK_zscore, age_at_vaccination){ 
  #tau * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl))

  #params
  TVCL <- 0.0051
  TVV1 <- 2.57
  TVV2 <- 0.944
  TVQ <- 0.172
  TVV3 <- 1.47
  TVQ2 <- 0.00782
  
  #Age discretisation
  age_v <- c(0.0,0.5,1,1.5,2,2.5,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18) * 365
  
  W_L_F <- c(0.0402000, -0.1507000, -0.2384000, -0.2815000, -0.3032000, -0.3155000, -0.3283000, -0.3440000,
             -0.4834000, -0.5185000, -0.5495000, -0.5740000, -0.5905000, -0.5735901, 
            -0.5587019, -0.6037183, -0.7121661, -0.8789837, -1.0778904, -1.2511913, -1.3269735, -1.3141719)
  W_M_F <- c(5.84580,  8.22540,  9.60080, 10.85340, 12.10150, 13.28370, 14.97270, 17.15510, 19.12760,
             21.22740, 23.63690, 26.55190, 29.96630,
              33.78409, 38.03438, 42.20747, 45.92315, 48.88839, 51.00794, 52.42775, 53.47799, 53.98051)
  W_S_F <- c(0.1261900, 0.1219900, 0.1229900, 0.1233500, 0.1247200, 0.1273700, 0.1337600, 0.1437100,
             0.1456900, 0.1523000, 0.1587600,
             0.1648300, 0.1702500, 0.1756898, 0.1785063, 0.1763689, 0.1701308, 0.1613282, 0.1520297, 0.1446439,
             0.1411725, 0.1411252)
  
  W_L_M <- c(0.1738000,  0.0917000,  0.0413000,  0.0029000, -0.0289000, -0.0564000, -0.0920000,
             -0.1325000, -0.2548000, -0.3804000,
              -0.4964000, -0.5946000, -0.6624000, -0.6358593, -0.5270620, -0.4333559, -0.3728613, -0.3619037,
             -0.4122345, -0.5075487, -0.5901793, -0.6077477)
  W_M_M <- c(6.37620,  8.90140, 10.31080, 11.54860, 12.74010, 13.83090, 15.34860, 17.34520, 19.39400,
             21.68100, 24.13710, 26.73580, 29.57360,
             32.98352, 37.14884, 41.89264, 47.05580, 52.29523, 57.13584, 61.14131, 64.12427, 65.27123)
  W_S_M <- c(0.1172700, 0.1088100, 0.1100700, 0.1126100, 0.1160400, 0.1195300, 0.1242500, 0.1313300,
             0.1317800, 0.1355400, 0.1401600, 0.1475200, 0.1576000, 0.1682108, 0.1765885, 0.1796924,
             0.1768534, 0.1696619, 0.1609673, 0.1534760, 0.1487292, 0.1473498)

  amt <- 700 # 

  vnapply(
    seq_along(tbv_iiva),
    function(i) {
      
      iiva <- tbv_iiva[[i]]
      iivb <- tbv_iivb[[i]]
      zscore <- tbv_PK_zscore[[i]]
      sex <- tbv_PK_sx[[i]]
      age_agent <- age_at_vaccination[[i]]
      tt <- t[[i]] #JDC: new addition
      
      age_index <- max(which(age_v <= age_agent))
      
      WT <- 0 # good way of indicating an issue?
      
      if(sex==1){
        WT <- W_M_M[age_index] *(1 + zscore*W_L_M[age_index]*W_S_M[age_index])^(1/W_L_M[age_index])
      }else{
        WT <- W_M_F[age_index] *(1 + zscore*W_L_F[age_index]*W_S_F[age_index])^(1/W_L_F[age_index])
      }
      
      WTCL <- (WT/70)**0.75
      WTV <- (WT/70)**1
      
      V2 <- TVV2*WTV
      Q <- TVQ*WTCL
      V3 <- TVV3*WTV
      Q2 <- TVQ2*WTCL
      
      CL <- TVCL*exp( iiva )*WTCL #IIV
      V1 <- TVV1*exp( iivb )*WTV #IIV
      #Rate consts.
      k10 <- 24*CL/V1
      k12 <- 24*Q/V1
      k21 <- 24*Q/V2
      k13 <- 24*Q2/V1
      k31 <- 24*Q2/V3
      
      mt <- matrix(c(-(k12+k10+k13), k21, k31, k12, -k21, 0, k13, 0, -k31), nrow = 3, byrow = T)
      p <- eigen(mt)$vectors #These are in the right order
      lambda <- eigen(mt)$values
      y0 <- c(amt,0,0) 
      zz0 <- solve(p)%*%y0
      z <- c(zz0[1]*exp(lambda[1] * tt ),
             zz0[2]*exp(lambda[2] * tt ),
             zz0[3]*exp(lambda[3] * tt) )
      zz <-  p[1,1]*z[1] + p[1,2]*z[2] + p[1,3]*z[3]# p%*%z
      zz/V1
    }
  )
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

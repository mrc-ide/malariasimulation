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
    antibodies <- calculate_tbv_antibodies(
      timestep - time_vaccinated[vaccinated_in_state]#,
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

#correlated r.v.s
chol_two <- function(v1,v2,c12){
  z1 <- rnorm(1,0,1)
  z2 <- rnorm(1,0,1)
  x1 <- sqrt(v1)*z1
  x2 <- c12*z1/sqrt(v1) + sqrt(-(c12**2)/v1 + v2)*z2
  return(c(x1,x2))
}

#RK4 functions
f1 <- function(x, y, z, K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31){
  K21r*y - K12r*x - K10r*x - K13r*x + K31r*z
}

f2 <- function(x, y, z, K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31){
  K12r*x - K21r*y
}

f3 <- function(x, y, z, K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31){
  K13r*x - K31r*z
}

calculate_tbv_antibodies <- function(t){ #this'll need updating. Add more arguments??
  #tau * (rho * exp(-t * log(2) / ds) + (1 - rho) * exp(-t * log(2) / dl))
  #199*(0.36*exp(-t/6) + 0.44*exp(-t/37) + 0.2*exp(-t/50))
  #Is there a problem if t near 0???
    #tt <- t#ifelse(t>0.3, t,0.3 )
    #  k10 <- (0.0051/2.57)*24
    #  k12 <- (0.172/2.57)*24
    #  k21 <- (0.172/0.944)*24
    #  k13 <- (0.00782/2.57)*24
    #  k31 <- (0.00782/1.47)*24
    #  amt <- 700/2.57
    #  mt <- matrix( c(-(k12+k10+k13)*t, k21*t, k31*t, k12*t, -k21*t, 0, k13*t, 0, -k31*t), nrow = 3,byrow=T) 
    #  me <- Matrix::expm(mt)
    # (me[1,1])*amt
   #
  #if(t > 350){ # might need increasing, due to IIV...
  #  0.0
  #} else{
  
  #tt <- t+1.7

  # #params
  # TVCL <- 0.0051
  # TVV1 <- 2.57
  # TVV2 <- 0.944
  # TVQ <- 0.172
  # TVV3 <- 1.47
  # TVQ2 <- 0.0782
  # WT <- 70 #fixed for us?
  # AMT <- WT*10
  # 
  # WTCL <- (WT/70)**0.75
  # WTV <- (WT/70)**1
  # 
  # #randos
  # #if(IIV==0){
  #   r1 <- 0
  #   r2 <- 0
  # #}else{
  # #  r <- chol_two(0.0864,0.117,0.0797)
  # #  r1 <- r[1]
  # #  r2 <- r[2]
  # #}
  # CL <- TVCL*exp(r1)*WTCL
  # V1 <- TVV1*exp(r2)*WTV
  # V2 <- TVV2*WTV
  # Q <- TVQ*WTCL
  # V3 <- TVV3*WTV
  # Q2 <- TVQ2*WTCL
  # 
  # #Rate consts.
  # K10 <- 24*CL/V1
  # K12 <- 24*Q/V1
  # K21 <- 24*Q/V2
  # K13 <- 24*Q2/V1
  # K31 <- 24*Q2/V3
  # 
  # k1 <- k2 <- k3 <- k4 <- rep(0,3)
  # 
  # #now a for loop. How long to run the model for? 'time' is in days
  # dt <- 0.2
  # TL <- max(as.integer((1/dt) * t ),1) # Think about this
  # #print(paste0('TL: ',TL))
  # 
  # #if(t<0.25){
  # #  AMT/V1
  # #}else{
  # 
  #   CENT <- rep(0, TL+1)
  #   PERI1 <- rep(0, TL+1)
  #   PERI2 <- rep(0, TL+1)
  # 
  #   CENT[1] <- AMT
  # 
  #   for(k in 1:TL){
  # 
  #     xL = CENT[k]
  #     yL = PERI1[k]
  #     zL = PERI2[k]
  # 
  #     k1[1] <- dt*f1(x=xL, y=yL, z=zL, K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  #     k1[2] <- dt*f2(x=xL, y=yL, z=zL, K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  #     k1[3] <- dt*f3(x=xL, y=yL, z=zL, K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  # 
  #     k2[1] <- dt*f1(x=xL + 0.5*k1[1], y=yL + 0.5*k1[2], z=zL + 0.5*k1[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  #     k2[2] <- dt*f2(x=xL + 0.5*k1[1], y=yL + 0.5*k1[2], z=zL + 0.5*k1[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  #     k2[3] <- dt*f3(x=xL + 0.5*k1[1], y=yL + 0.5*k1[2], z=zL + 0.5*k1[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  # 
  #     k3[1] <- dt*f1(x=xL + 0.5*k2[1], y=yL + 0.5*k2[2], z=zL + 0.5*k2[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  #     k3[2] <- dt*f2(x=xL + 0.5*k2[1], y=yL + 0.5*k2[2], z=zL + 0.5*k2[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  #     k3[3] <- dt*f3(x=xL + 0.5*k2[1], y=yL + 0.5*k2[2], z=zL + 0.5*k2[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  # 
  #     k4[1] <- dt*f1(x=xL + k3[1], y=yL + k3[2], z=zL + k3[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  #     k4[2] <- dt*f2(x=xL + k3[1], y=yL + k3[2], z=zL + k3[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  #     k4[3] <- dt*f3(x=xL + k3[1], y=yL + k3[2], z=zL + k3[3], K10r=K10, K12r=K12, K21r=K21, K13r=K13, K31r=K31)
  # 
  #     CENT[k+1] = CENT[k] +  (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])/6
  #     PERI1[k+1] = PERI1[k] +  (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2])/6
  #     PERI2[k+1] = PERI2[k] +  (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3])/6
  # 
  #   } #time loop
  # 
  #   #conc <-
  #   CENT[TL+1]/V1
  #}
  #} #what were these two for?
  
  ##Now diag
  #params
  TVCL <- 0.0051
  TVV1 <- 2.57
  TVV2 <- 0.944
  TVQ <- 0.172
  TVV3 <- 1.47
  TVQ2 <- 0.0782
  WT <- 70 #fixed for us?
  AMT <- WT*10

  WTCL <- (WT/70)**0.75
  WTV <- (WT/70)**1

  #randos
  #if(IIV==0){
    r1 <- 0
    r2 <- 0
  #}else{
  #  r <- chol_two(0.0864,0.117,0.0797)
  #  r1 <- r[1]
  #  r2 <- r[2]
  #}
  CL <- TVCL*exp(r1)*WTCL
  V1 <- TVV1*exp(r2)*WTV
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

  amt <- 700 # mg/kg * WT

  mt <- matrix(c(-(k12+k10+k13), k21, k31, k12, -k21, 0, k13, 0, -k31), nrow = 3, byrow = T)
  p <- eigen(mt)$vectors #These in the right order
  lambda <- eigen(mt)$values
  y0 <- c(amt,0,0) 
  zz0 <- solve(p)%*%y0
  z <- c(zz0[1]*exp(lambda[1] * max( t ,.1)), 
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

calculate_TBA <- function(mx, k, tra) { #same
  offset <- (k / (k + mx)) ^ k
  scale <- 1. / (1. - offset)
  tra_transformation <- (k / (k + mx * (1 - tra))) ^ k
  scale * (tra_transformation - offset)
}
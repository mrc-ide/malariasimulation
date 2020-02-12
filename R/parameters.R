#' @title Get model parameters
#' @description
#' get_paramaters creates a named list of parameters for use in the model. These
#' parameters are passed to process functions. These parameters are explained in
#' "The US President's Malaria Initiative, Plasmodium falciparum transmission
#' and mortality: A modelling study."
#'
#' NOTE: this function is likely to be extended to read in command line / config
#' file parameters
#'
#' The parameters are defined below.
#'
#' fixed state transitions:
#'
#' * rd - the rate at which humans move from state D to A
#' * ra - the rate at which humans move from state A to U
#' * ru - the rate at which humans move from state U to S
#' * rel - the rate at which mosquitos move from state E to L
#' * rl - the rate at which mosquitos move from state L to P
#' * rpl - the rate at which mosquitos move from state P to Sm
#' * mup - the rate at which pupal mosquitos die
#' * mum - the rate at which developed mosquitos die
#'
#' immunity decay rates:
#'
#' * rm - decay rate for maternal immunity to clinical disease
#' * rvm - decay rate for maternal immunity to severe disease
#' * rb - decay rate for acquired pre-erytrhrocytic immunity
#' * rc - decay rate for acquired immunity to clinical disease
#' * rva - decay rate for acquired immunity to severe disease
#' * rid - decay rate for acquired immunity to detectability
#'
#' probability of pre-erythrocytic infection:
#'
#' * b0 - maximum probability due to no immunity
#' * b1 - maximum reduction due to immunity
#' * ib0 - scale parameter
#' * kb - shape parameter
#'
#' probability of clinical infection:
#'
#' * phi0 - maximum probability due to no immunity
#' * phi1 - maximum reduction due to immunity
#' * ic0 - scale parameter
#' * kc - shape parameter
#'
#' probability of severe infection:
#'
#' * theta0 - maximum probability due to no immunity
#' * theta1 - maximum reduction due to immunity
#' * iv0 - scale parameter
#' * kv - shape parameter
#' * fv0 - age dependent modifier
#' * av - age dependent modifier
#' * gammav - age dependent modifier
#'
#' immunity reducing probability of detection:
#' * fd0 - time-scale at which immunity changes with age
#' * ad - scale parameter relating age to immunity
#' * gammad - shape parameter relating age to immunity
#' * d1 - minimum probability due to immunity
#' * dmin - minimum probability due to immunity NOTE: there appears to be a
#' mistake here!
#' * id0 - scale parameter 
#' * kd - shape parameter
#'
#' immunity boost grace periods:
#'
#' * ub - period in which pre-erythrocytic immunity cannot be boosted
#' * uc - period in which clinical immunity cannot be boosted
#' * uv - period in which severe immunity cannot be boosted
#' * ud - period in which immunity to detectability cannot be boosted
#'
#' blood meal rates:
#'
#' * av1 - blood meal rate for the first variety of mosquitos
#' * av2 - blood meal rate for the second variety of mosquitos
#' * av3 - blood meal rate for the third variety of mosquitos
#'
#' infectivity towards mosquitos:
#'
#' * cd - infectivity of clinically diseased humans towards mosquitos
#' * gamma1- parameter for infectivity of asymptomatic humans
#' * cu - infectivity of sub-patent infection
#'
#' unique biting rate:
#'
#' * a0 - age dependent biting parameter
#' * rho - age dependent biting parameter
#' * sigma_squared - heterogeneity parameter
#'
#' miscellaneous:
#'
#' * de - delay for infection
#' * beta - the average number of eggs laid per female mosquito per day
#' * human population - the number of humans to model
#' * mosquito limit - the maximum number of mosquitos to allow for in the
#' * days_per_timestep - the number of days to model per timestep
get_parameters <- function() {
  days_per_timestep <- 1
  human_population <- 100 * 1000
  parameters <- list(
    rd    = days_per_timestep / 5,
    ra    = days_per_timestep / 195,
    ru    = days_per_timestep / 110,
    rel   = days_per_timestep / 6.64,
    rl    = days_per_timestep / 3.72,
    rpl   = days_per_timestep / .643,
    mup   = days_per_timestep * .249,
    mum   = days_per_timestep * .249, #NOTE: set from sitefile
    sigma_squared   = 1.67,
    n_heterogeneity_groups = 5,
    # immunity decay rates
    rm    = days_per_timestep / 67.6952,
    rvm   = days_per_timestep / 76.8365,
    rb    = days_per_timestep / (10 * 365),
    rc    = days_per_timestep / (30 * 365),
    rva   = days_per_timestep / (30 * 365),
    rid   = days_per_timestep / (10 * 365),
    # blood immunity parameters
    b0    = 0.590076,
    b1    = 0.5,
    ib0   = 43.8787,
    kb    = 2.15506,
    # immunity boost grace periods
    ub    = 7.19919,
    uc    = 67.6952,
    uv    = 11.4321,
    ud    = 9.44512,
    # blood meal rates
    av1   = .92,
    av2   = .74,
    av3   = .94,
    # infectivity towards mosquitos
    cd    = 0.068,
    gamma1= 1.82425,
    cu    = 0.00062,
    # unique biting rate
    a0    = 8 * 365 / days_per_timestep,
    rho   = .85,
    # clinical immunity parameters
    phi0  = .0749886,
    phi1  = .0001191,
    ic0   = 18.02366,
    kc    = 2.36949,
    # severe disease immunity parameters
    theta0  = .0749886,
    theta1  = .0001191,
    kv      = 2.00048,
    fv0     = 0.141195,
    av      = 2493.41,
    gammav  = 2.91282,
    iv0     = 1.09629,
    # delay for infection
    de      = 12 / days_per_timestep,
    # asymptomatic immunity parameters
    fd0   = 0.007055,
    ad    = 21.9 * 365 / days_per_timestep,
    gammad= 4.8183,
    d1    = 0.160527,
    dmin  = 0, #NOTE: what should this be?
    id0   = 1.577533,
    kd    = .476614,
    # mortality parameters
    ftv   = .5,
    mortality_probability_table = rep(.05, 100),
    v     = .065, # NOTE: there are two definitions of this: one on line 124 and one in the parameters table
    pcm   = .774368,
    pvm   = .195768,
    # carrying capacity parameters
    K0    = 10,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    gamma = 13.25,
    # larval mortality rates
    me    = days_per_timestep * .0338,
    ml    = days_per_timestep * .0348,
    # egg laying parameter
    beta  = days_per_timestep * 21.2,
    human_population = human_population,
    mosquito_limit   = 100 * human_population,
    days_per_timestep  = days_per_timestep
  )

	parameters$R_bar <- mean(vnapply(1:365, function(t) rainfall(
		t,
		parameters$days_per_timestep,
    parameters$g0,
		c(parameters$g1, parameters$g2, parameters$g3),
		c(parameters$h1, parameters$h2, parameters$h3)
	)))

  parameters
}

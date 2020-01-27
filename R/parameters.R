#' @description
#'
#' get_paramaters creates a list of parameters for use in the model. These
#' parameters are passed to process functions
#'
#' NOTE: this function is likely to be extended to read in command line / config
#' file parameters
get_parameters <- function() {
  timestep_to_day <- 1
  parameters <- list(
    b0    = 0.590076,
    b1    = 0.5,
    ib0   = 43.8787,
    kb    = 2.15506,
    rd    = 1/5,
    ra    = 1/195,
    ru    = 1/110,
    rt    = 1/5,
    ft    = 1/2, # NOTE: set from sitefile
    av1   = .92,
    av2   = .74,
    av3   = .94,
    cd    = 0.068,
    ct    = 0.021896,
    gamma1= 1.82425,
    cu    = 0.00062,
    rm    = 1 / (67.6952 * timestep_to_day),
    rvm   = 1 / (76.8365 * timestep_to_day),
    rb    = 1 / (10 * 365 * timestep_to_day),
    rc    = 1 / (30 * 365 * timestep_to_day),
    rva   = 1 / (30 * 365 * timestep_to_day),
    rd    = 1 / (10 * 365 * timestep_to_day),
    ub    = 7.19919,
    uc    = 67.6952,
    uv    = 11.4321,
    ud    = 9.44512,
    a0    = 8 * 365 * timestep_to_day,
    phi0  = .0749886,
    phi1  = .0001191,
    ic0   = 18.02366,
    kc    = 2.36949,
    rho   = .85,
    sigma_squared   = 1.67,
    theta0  = .0749886,
    theta1  = .0001191,
    kv      = 2.00048,
    fv0     = 0.141195,
    av      = 2493.41,
    gammav  = 2.91282,
    iv0     = 1.09629,
    de      = 12 * timestep_to_day,
    fd0   = 0.007055,
    ad    = 21.9 * 365 * timestep_to_day,
    gammad= 4.8183,
    d1    = 0.160527,
    dmin  = 0, #NOTE: what should this be?
    id0   = 1.577533,
    kd    = .476614,
    ftv   = .5,
    pcm   = .774368,
    pvm   = .195768,
    v     = .065, # NOTE: there are two definitions of this: one on line 124 and one in the parameters table
    rel   = 1 / (6.64 * timestep_to_day),
    rl    = 1 / (3.72 * timestep_to_day),
    rpl   = 1 / (.643 * timestep_to_day),
    beta  = 21.2,
    K0    = 10,
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    mup   = .249,
    mum   = .249, #NOTE: set from sitefile
    me    = .0338,
    ml    = .0348,
    gamma = 13.25,
    mortality_probability_table = rep(.05, 100),
    timestep_to_day = timestep_to_day
  )

	parameters$R_bar <- mean(vnapply(1:365, function(t) rainfall(
		t,
		parameters$timestep_to_day,
    parameters$g0,
		c(parameters$g1, parameters$g2, parameters$g3),
		c(parameters$h1, parameters$h2, parameters$h3)
	)))

  parameters
}

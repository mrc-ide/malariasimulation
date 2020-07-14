#' @title Get model parameters
#' @description
#' get_paramaters creates a named list of parameters for use in the model. These
#' parameters are passed to process functions. These parameters are explained in
#' "The US President's Malaria Initiative, Plasmodium falciparum transmission
#' and mortality: A modelling study."
#'
#' @param overrides a named list of parameter values to use instead of defaults
#'
#' NOTE: this function is likely to be extended to read in command line / config
#' file parameters
#'
#' The parameters are defined below.
#'
#' fixed state transitions:
#'
#' * dd - the delay for humans to move from state D to A
#' * da - the delay for humans to move from state A to U
#' * du - the delay for humans to move from state U to S
#' * del - the delay for mosquitos to move from state E to L
#' * dl - the delay for mosquitos to move from state L to P
#' * dpl - the delay mosquitos to move from state P to Sm
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
#' * severe_enabled - whether to model severe disease
#' * theta0 - maximum probability due to no immunity
#' * theta1 - maximum reduction due to immunity
#' * iv0 - scale parameter
#' * kv - shape parameter
#' * fv0 - age dependent modifier
#' * av - age dependent modifier
#' * gammav - age dependent modifier
#'
#' immunity reducing probability of detection:
#'
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
#' infectivity towards mosquitos:
#'
#' * cd - infectivity of clinically diseased humans towards mosquitos
#' * gamma1 - parameter for infectivity of asymptomatic humans
#' * cu - infectivity of sub-patent infection
#'
#' unique biting rate:
#'
#' * a0 - age dependent biting parameter
#' * rho - age dependent biting parameter
#' * sigma_squared - heterogeneity parameter
#' * n_heterogeneity_groups - number discretised groups for heterogeneity, used
#' for sampling mothers
#'
#' mortality parameters:
#'
#' * mortality_rate - human mortality rate across age groups
#' * v - mortality scaling factor from severe disease
#' * pcm - new-born clinical immunity relative to mother's
#' * pvm - new-born severe immunity relative to mother's
#' * me - early stage larval mortality rate
#' * ml - late stage larval mortality rate
#'
#' carrying capacity parameters:
#'
#' * model_seasonality - boolean switch TRUE iff the simulation models seasonal rainfall
#' * K0 - carrying capacity (derived)
#' * g0 to g3 - rainfall shape parameters
#' * h1 to h3 - rainfall shape parameters
#' * gamma - effect of density dependence on late instars relative to early
#' instars
#'
#' initial state proportions:
#'
#' * s_proportion - the proportion of `human_population` that begin as Susceptable
#' * d_proportion - the proportion of `human_population` that begin with
#' clinical disease
#' * a_proportion - the proportion of `human_population` that begin as
#' Asymptomatic
#' * u_proportion - the proportion of `human_population` that begin as
#' Subpatents
#'
#' initial immunity values:
#'
#' * init_icm - the immunity from clinical disease at birth
#' * init_ivm - the immunity from severe disease at birth
#' * init_ib  - the initial pre-erythrocitic immunity
#' * init_ica - the initial acquired immunity from clinical disease
#' * init_iva - the initial acquired immunity from severe disease
#' * init_id  - the initial acquired immunity to detectability
#'
#' incubation periods:
#'
#' * de - delay for infection
#' * dem - delay for infection in mosquitoes
#'
#' vector biology:
#'
#' * total_M - the initial number of adult mosquitoes
#' * beta - the average number of eggs laid per female mosquito per day
#' * variety_proportions - the relative proportions of each species
#' * blood_meal_rates - the blood meal rates for each species
#'
#' miscellaneous:
#'
#' * human_population - the number of humans to model
#' * mosquito_limit - the maximum number of mosquitos to allow for in the
#' simulation
#' * days_per_timestep - the number of days to model per timestep
#' * vector_ode - whether to use the ODE model to model mosquitos
#' @export
get_parameters <- function(overrides = list()) {
  days_per_timestep <- 1
  human_population <- 100
  if ('human_population' %in% names(overrides)) {
    human_population <- overrides$human_population
  }
  parameters <- list(
    dd    = 5,
    da    = 195,
    du    = 110,
    del   = 6.64,
    dl    = 3.72,
    dpl   = .643,
    mup   = .249,
    mum   = .249, #NOTE: set from sitefile
    sigma_squared   = 1.67,
    n_heterogeneity_groups = 5,
    # immunity decay rates
    rm    = 67.6952,
    rvm   = 76.8365,
    rb    = 10 * 365,
    rc    = 30 * 365,
    rva   = 30 * 365,
    rid   = 10 * 365,
    # blood immunity parameters
    b0    = 0.590076,
    b1    = 0.5,
    ib0   = 43.8787,
    kb    = 2.15506,
    # immunity boost grace periods
    ub    = 7.2,
    uc    = 6.06,
    uv    = 11.4321,
    ud    = 9.44512,
    # infectivity towards mosquitos
    cd    = 0.068,
    gamma1= 1.82425,
    cu    = 0.0062,
    # unique biting rate
    a0    = 8 * 365,
    rho   = .85,
    # clinical immunity parameters
    phi0  = .0749886,
    phi1  = .0001191,
    ic0   = 18.02366,
    kc    = 2.36949,
    # severe disease immunity parameters
    severe_enabled = 0,
    theta0  = .0749886,
    theta1  = .0001191,
    kv      = 2.00048,
    fv0     = 0.141195,
    av      = 2493.41,
    gammav  = 2.91282,
    iv0     = 1.09629,
    # delay for infection
    de      = 12,
    dem     = 10,
    # asymptomatic immunity parameters
    fd0   = 0.007055,
    ad    = 21.9 * 365,
    gammad= 4.8183,
    d1    = 0.160527,
    id0   = 1.577533,
    kd    = .476614,
    # mortality parameters
    average_age = 7663 / days_per_timestep,
    v     = .065, # NOTE: there are two definitions of this: one on line 124 and one in the parameters table
    pcm   = .774368,
    pvm   = .195768,
    # carrying capacity parameters
    g0    = 2,
    g1   = .3,
    g2   = .6,
    g3   = .9,
    h1   = .1,
    h2   = .4,
    h3   = .7,
    gamma = 13.25,
    model_seasonality = FALSE,
    # larval mortality rates
    me    = .0338,
    ml    = .0348,
    # initial state proportions
    s_proportion = 0.420433246,
    d_proportion = 0.007215064,
    a_proportion = 0.439323667,
    u_proportion = 0.133028023,
    # initial immunities
    init_ica = 0,
    init_iva = 0,
    init_icm = 0,
    init_ivm = 0,
    init_id  = 0,
    init_ib  = 0,
    # vector biology
    beta     = 21.2,
    total_M  = 1000,
    variety_proportions = c(.5, .3, .2),
    blood_meal_rates    = c(.92, .74, .94),
    # misc
    human_population = human_population,
    mosquito_limit   = 10000 * human_population,
    vector_ode         = FALSE,
    days_per_timestep  = days_per_timestep
  )

  parameters$mortality_rate = 1 - exp(
    -days_per_timestep * (1 / parameters$average_age)
  )

	parameters$R_bar <- mean(vnapply(1:365, function(t) rainfall(
		t,
		parameters$days_per_timestep,
    parameters$g0,
		c(parameters$g1, parameters$g2, parameters$g3),
		c(parameters$h1, parameters$h2, parameters$h3)
	)))

  # Override parameters with any client specified ones
  if (!is.list(overrides)) {
    stop('overrides must be a list')
  }

  for (name in names(overrides)) {
    if (!(name %in% names(parameters))) {
      stop(paste('unknown parameter', name, sep=' '))
    }
    parameters[[name]] <- overrides[[name]]
  }

  props <- c(
    parameters$s_proportion,
    parameters$d_proportion,
    parameters$a_proportion,
    parameters$u_proportion
  )

  if (!all.equal(sum(props), 1)) {
    stop("Starting proportions do not sum to 1")
  }

  parameters$K0 <- calculate_carrying_capacity(parameters)

  parameters
}

#' @title Parameterise equilibrium proportions
#' @description parameterise equilibrium proportions from a list
#'
#' @param parameters a named list of model parameters
#' @param state_props the equilibrium proportions in a named list
#' @export
parameterise_equilibrium <- function(parameters, state_props) {
  parameters$s_proportion <- state_props[['S']]
  parameters$d_proportion <- state_props[['D']]
  parameters$a_proportion <- state_props[['A']]
  parameters$u_proportion <- state_props[['U']]
  parameters
}

#' @title Get model parameters
#' @description
#' get_paramaters creates a named list of parameters for use in the model. These
#' parameters are passed to process functions. These parameters are explained in
#' "The US President's Malaria Initiative, Plasmodium falciparum transmission
#' and mortality: A modelling study."
#'
#' @param overrides a named list of parameter values to use instead of defaults
#' The parameters are defined below.
#'
#' fixed state transitions:
#'
#' * dd - the delay for humans to move from state D to A
#' * dt - the delay for humans to move from state Tr to Ph
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
#' * fvt - reduced probability of death due to treatment
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
#' * ct - infectivity of treated infection
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
#' * average_age - the average age of humans (in timesteps)
#' * v - mortality scaling factor from severe disease
#' * pcm - new-born clinical immunity relative to mother's
#' * pvm - new-born severe immunity relative to mother's
#' * me - early stage larval mortality rate
#' * ml - late stage larval mortality rate
#'
#' carrying capacity parameters:
#'
#' * model_seasonality - boolean switch TRUE iff the simulation models seasonal rainfall
#' * g0 - rainfall fourier parameter
#' * g - rainfall fourier parameter
#' * h - rainfall fourier parameters
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
#' subpatents
#' * t_proportion - the proportion of `human_population` that begin treated
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
#' species specific values are vectors
#'
#' * beta - the average number of eggs laid per female mosquito per day
#' * total_M - the initial number of adult mosquitos in the simulation
#' * init_foim - the FOIM used to calculate the equilibrium state for mosquitoes
#' * species - names of the species in the simulation
#' * species_proportions - the relative proportions of each species
#' * blood_meal_rates - the blood meal rates for each species
#' * Q0 - proportion of blood meals taken on humans
#' * foraging_time - time spent taking blood meals
#'
#' feeding cycle:
#'
#' * bednets - boolean for if bednets are enabled
#' * rn - probability mosquito is repelled by the bednet
#' * rnm - minimum probability mosquito is repelled by the bednet 
#' * dn0 - probability killed by the bednet
#' * spraying - boolean for if indoor spraying is enabled
#' * rs - probability repelled by indoor spraying
#' * phi_indoors - proportion of bites taken indoors
#' * phi_bednets - proportion of bites taken in bed
#' * endophily - proportion of mosquitoes resting indoors after feeding with no
#' intervention
#' * gammas - the half-life of spraying efficacy (timesteps)
#' * gamman - the half-life of bednet efficacy (timesteps)
#'
#' please set vector control strategies using `set_betnets` and `set_spraying`
#'
#' treatment parameters:
#' I recommend setting these with the convenience functions in
#' `drug_parameters.R`
#'
#' * drug_efficacy - a vector of efficacies for available drugs
#' * drug_rel_c - a vector of relative onwards infectiousness values for drugs
#' * drug_prophylaxis_shape - a vector of shape parameters for weibull curves to
#' model prophylaxis for each drug
#' * drug_prophylaxis_scale - a vector of scale parameters for weibull curves to
#' model prophylaxis for each drug
#' * ft - probability of seeking treatment if clinically diseased
#' * clinical_treatment_drugs - a vector of drugs that are avaliable for
#' clinically diseased (these values refer to the index in drug_* parameters)
#' * clinical_treatment_coverage - a vector of coverage values for each drug
#'
#' RTS,S paramters:
#'
#' * rtss_vmax - the maximum efficacy of the vaccine
#' * rtss_alpha - shape parameter for the vaccine efficacy model
#' * rtss_beta - scale parameter for the vaccine efficacy model
#' * rtss_cs - peak parameters for the antibody model (mean and std. dev)
#' * rtss_cs_boost - peak parameters for the antibody model for booster rounds (mean and std. dev)
#' * rtss_rho - delay parameters for the antibody model (mean and std. dev)
#' * rtss_rho_boost - delay parameters for the antibody model for booster rounds (mean and std. dev)
#' * rtss_ds - delay parameters for the antibody model (mean and std. dev)
#' * rtss_dl - delay parameters for the antibody model (mean and std. dev)
#'
#' I recommend setting strategies with the convenience functions in
#' `vaccine_parameters.R:set_rtss`
#'
#' MDA and SMC parameters:
#' We recommend setting these with convenience functions in `mda_parameters.R`
#'
#' TBV parameters:
#'
#' * tbv_mt - effect on treated infectiousness
#' * tbv_md - effect on diseased infectiousness
#' * tbv_ma - effect on asymptomatic infectiousness
#' * tbv_mu - effect on subpatent infectiousness
#' * tbv_k  - scale parameter for effect on infectiousness
#' * tbv_tau - peak antibody parameter
#' * tbv_rho - antibody component parameter
#' * tbv_ds - antibody short-term delay parameter
#' * tbv_dl - antibody long-term delay parameter
#' * tbv_tra_mu - transmission reduction parameter
#' * tbv_gamma1 - transmission reduction parameter
#' * tbv_gamma2 - transmission reduction parameter
#'
#' I recommend setting tbv strategies with the convenience functions in
#' `vaccine_parameters.R`, these are the same as for RTS,S
#'
#' rendering:
#' All values are in timesteps and all ranges are inclusive
#'
#' * prevalence_rendering_min_ages - the minimum ages for clinical prevalence
#' outputs
#' * prevalence_rendering_max_ages - the corresponding max ages
#' * incidence_rendering_min_ages - the minimum ages for clinical incidence
#' outputs
#' * incidence_rendering_max_ages - the corresponding max ages
#' * severe_prevalence_rendering_min_ages - the minimum ages for severe
#' prevalence outputs
#' * severe_prevalence_rendering_max_ages - the corresponding max ages
#' * severe_incidence_rendering_min_ages - the minimum ages for severe incidence
#' outputs
#' * severe_incidence_rendering_max_ages - the corresponding max ages
#'
#' miscellaneous:
#'
#' * human_population - the number of humans to model
#' * mosquito_limit - the maximum number of mosquitos to allow for in the
#' simulation
#' * days_per_timestep - the number of days to model per timestep
#'
#' @export
get_parameters <- function(overrides = list()) {
  days_per_timestep <- 1
  parameters <- list(
    dd    = 5,
    dt    = 5,
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
    ct    = 0.021896,
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
    fvt     = 0.5,
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
    v     = .065, # NOTE: there are two definitions of this in the literature: one on line 124 and one in the parameters table
    pcm   = .774368,
    pvm   = .195768,
    # carrying capacity parameters
    g0    = 2,
    g     = c(.3, .6, .9),
    h     = c(.1, .4, .7),
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
    t_proportion = 0,
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
    init_foim= 0,
    # order of species: An gambiae s.s, An arabiensis, An funestus
    species             = c('gamb', 'arab', 'fun'),
    species_proportions = c(.5, .3, .2),
    blood_meal_rates    = rep(1/3, 3),
    Q0                  = c(.92, .71, .94),
    endophily           = c(.813, .422, .813),
    foraging_time       = .69,
    # bed nets
    bednets = FALSE,
    rn = c(.56, .46, .56),
    rnm = c(.24, .1, .24),
    dn0 = rep(.533, 3),
    phi_bednets = c(.89, .9, .9),
    gamman = 2.64 * 365,
    # indoor spraying
    spraying = FALSE,
    rs = rep(.2, 3),
    phi_indoors = c(.97, .96, .98),
    gammas = .25 * 365,
    # treatment
    drug_efficacy          = numeric(0),
    drug_rel_c             = numeric(0),
    drug_prophilaxis_shape = numeric(0),
    drug_prophilaxis_scale = numeric(0),
    clinical_treatment_drugs     = numeric(0),
    clinical_treatment_coverages = numeric(0),
    ft = 0,
    # rts,s
    rtss = FALSE,
    rtss_vmax = .93,
    rtss_alpha = .74,
    rtss_beta = 99.4,
    rtss_cs = c(6.37008, 0.35),
    rtss_cs_boost = c(5.56277, 0.35),
    rtss_rho = c(2.37832, 1.00813),
    rtss_rho_boost = c(1.03431, 1.02735),
    rtss_ds = c(3.74502, 0.341185),
    rtss_dl = c(6.30365, 0.396515),
    rtss_timeteps = NULL,
    rtss_coverages = NULL,
    rtss_ages = NULL,
    # MDA
    mda = FALSE,
    mda_drug = 0,
    mda_timesteps = NULL,
    mda_coverages = NULL,
    mda_min_age = -1,
    mda_max_age = -1,
    smc = FALSE,
    smc_drug = 0,
    smc_timesteps = NULL,
    smc_coverages = NULL,
    smc_min_age = -1,
    smc_max_age = -1,
    # tbv
    tbv = FALSE,
    tbv_mt = 35,
    tbv_md = 46.7,
    tbv_ma = 3.6,
    tbv_mu = 0.8,
    tbv_k = 0.9,
    tbv_tau = 22,
    tbv_rho = .7,
    tbv_ds = 45,
    tbv_dl = 591,
    tbv_tra_mu = 12.63,
    tbv_gamma1 = 2.5,
    tbv_gamma2 = .06,
    tbv_timeteps = NULL,
    tbv_coverages = NULL,
    tbv_frequency = -1,
    tbv_ages = NULL,
    # rendering
    prevalence_rendering_min_ages = 2 * 365,
    prevalence_rendering_max_ages = 10 * 365,
    incidence_rendering_min_ages = numeric(0),
    incidence_rendering_max_ages = numeric(0),
    severe_prevalence_rendering_min_ages = numeric(0),
    severe_prevalence_rendering_max_ages = numeric(0),
    severe_incidence_rendering_min_ages = numeric(0),
    severe_incidence_rendering_max_ages = numeric(0),
    # misc
    human_population = 100,
    mosquito_limit   = 100 * 1000,
    days_per_timestep  = days_per_timestep
  )

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
    parameters$u_proportion,
    parameters$t_proportion
  )

  if (!approx_sum(props, 1)) {
    stop("Starting proportions do not sum to 1")
  }

  parameters
}

#' @title Parameterise equilibrium proportions
#' @description parameterise equilibrium proportions from a list
#'
#' @param parameters the model parameters
#' @param eq the equilibrium solution output
#' @export
parameterise_human_equilibrium <- function(parameters, eq) {
  state_props <- colSums(eq$states[,c('S', 'D', 'A', 'U')])
  parameters$s_proportion <- state_props[['S']]
  parameters$d_proportion <- state_props[['D']]
  parameters$a_proportion <- state_props[['A']]
  parameters$u_proportion <- state_props[['U']]
  parameters$init_ica <- eq$states[,'ICA']
  parameters$init_icm <- eq$states[,'ICM']
  parameters$init_ib <- eq$states[,'IB']
  parameters$init_id <- eq$states[,'ID']
  parameters$init_foim <- eq$FOIM
  parameters
}

#' @title Parameterise total_M and carrying capacity for mosquitos from EIR
#'
#' @description NOTE: the inital EIR is likely to change unless the rest of the
#' model is in equilibrium. NOTE: please set seasonality first, since the mosquito_limit
#' will estimate an upper bound from the peak season.
#'
#' max_total_M is calculated using the equilibrium solution from "Modelling the
#' impact of vector control interventions on Anopheles gambiae population
#' dynamics"
#'
#' @param parameters the parameters to modify
#' @param EIR to work from
#' @export
parameterise_mosquito_equilibrium <- function(parameters, EIR) {
  parameterise_total_M(parameters. equilibrium_total_M(parameters, EIR))
}

#' @title Parameterise total_M
#'
#' @description Sets total_M and an upper bound for the number of mosquitoes in
#' the simulation. NOTE: please set seasonality first, since the mosquito_limit
#' will estimate an upper bound from the peak season.
#'
#' @param parameters the parameters to modify
#' @param total_M the initial adult mosquitoes in the simulation
#' @export
parameterise_total_M <- function(parameters, total_M) {
  parameters$total_M <- total_M
  K0 <- calculate_carrying_capacity(parameters)
  R_bar <- calculate_R_bar(parameters)
  max_K <- max(vnapply(seq(365), function(t) {
    carrying_capacity(
      t,
      parameters$model_seasonality,
      parameters$days_per_timestep,
      parameters$g0,
      parameters$g,
      parameters$h,
      K0,
      R_bar
    )
  }))
  omega <- calculate_omega(parameters)
  max_total_M <- max_K * (
    1 / (
      2 * parameters$dl * parameters$mum * (
        1 + parameters$dpl * parameters$mup
      )
    )
  ) * (
    1 / (
      parameters$gamma * (omega + 1)
    )
  ) * (
    omega / (parameters$ml * parameters$del) - (1 / (parameters$ml * parameters$dl)) - 1
  )
  parameters$mosquito_limit <- ceiling(max_total_M * 5) #Allow for random fluctuations
  parameters
}

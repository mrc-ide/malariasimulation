#' @title Get model parameters
#' @description
#' get_parameters creates a named list of parameters for use in the model. These
#' parameters are passed to process functions. The default parameters, which
#' are for Plasmodium falciparum, are explained in
#' "The US President's Malaria Initiative, Plasmodium falciparum transmission
#' and mortality: A modelling study."
#' 
#' Plasmodium vivax specific parameters are explained in
#' "Mathematical modelling of the impact of expanding levels of malaria control
#' interventions on Plasmodium vivax." by White, Michael T., et al.
#'
#' @param overrides a named list of parameter values to use instead of defaults
#' @param parasite Plasmodium parasite species ("falciparum" or "vivax"); default = "falciparum"
#' 
#' The parameters are defined below.
#'
#' human fixed state transitions:
#'
#' * dd - the delay for humans to move from state D to A; default = 5
#' * dt - the delay for humans to move from state Tr to Ph; default = 5
#' * da - the delay for humans to move from state A to U; default = 195
#' * du - the delay for humans to move from state U to S (P.f only); default = 110
#'
#' duration of subpatent infection: du (P.v only):
#' 
#' * du_max - Maximum duration of subpatent infection: default = 70
#' * du_min - Minimum duration of subpatent infection: default = 10
#' * ku - Shape parameter: default = 4.602
#' * au50 - Scale parameter: default = 9.9
#'
#' initial immunity values:
#'
#' * init_ica - the initial acquired immunity from clinical disease; default = 0
#' * init_icm - the immunity from clinical disease at birth; default = 0
#' * init_iva - the initial acquired immunity from severe disease (P.f only); default = 0
#' * init_ivm - the immunity from severe disease at birth (P.f only); default = 0
#' * init_ib  - the initial pre-erythrocitic immunity (P.f only); default = 0
#' * init_id  - the initial acquired immunity to detectability; default = 0
#' 
#' immunity boost grace periods:
#'
#' * ub - period in which pre-erythrocytic immunity cannot be boosted (P.f only); default = 7.2
#' * ud - period in which immunity to detectability cannot be boosted; default = 9.44512
#' * uc - period in which clinical immunity cannot be boosted; default = 6.06
#' * uv - period in which severe immunity cannot be boosted (P.f only); default = 11.4321
#'
#' maternal immunity parameters:
#' 
#' * pcm - new-born clinical immunity relative to mother's; default = 0.774368
#' * pvm - new-born severe immunity relative to mother's (P.f only); default = 0.195768
#' 
#' immunity decay rates:
#'
#' * rm - decay rate for maternal immunity to clinical disease; default = 67.6952
#' * rvm - decay rate for maternal immunity to severe disease (P.f only); default = 76.8365
#' * rb - decay rate for acquired pre-erythrocytic immunity (P.f only); default = 3650
#' * rid - decay rate for acquired immunity to detectability; default = 3650
#' * rc - decay rate for acquired immunity to clinical disease; default = 10950
#' * rva - decay rate for acquired immunity to severe disease (P.f only); default = 10950
#'
#' probability of pre-erythrocytic infection/blood immunity:
#'
#' * b0 - maximum probability due to no immunity (P.f only); default = 0.59
#' * b1 - maximum reduction due to immunity (P.f only); default = 0.5
#' * ib0 - scale parameter (P.f only); default = 43.9
#' * kb - shape parameter (P.f only); default = 2.16
#' * b - probability of pre-erythrocytic infection (P.v only): default = 0.5
#'
#' probability of asymptomatic detection/detectbility immunity (P.f only): 
#'
#' * fd0 - time-scale at which immunity changes with age; default = 0.007055
#' * ad - scale parameter relating age to immunity; default = 7993.5
#' * gammad - shape parameter relating age to immunity; default = 4.8183
#' * d1 - minimum probability due to immunity; default = 0.160527
#' * id0 - scale parameter; default = 1.577533
#' * kd - shape parameter; default = 0.476614
#'
#' probability of asymptomatic infection/antiparasite immunity (P.v only):
#'
#' * phi0lm - maximum probability due to no immunity; default = 0.8918
#' * phi1lm - maximum reduction due to immunity; default = 0.00482170890334156
#' * ic0lm - scale parameter; default = 27.52
#' * kclm - shape parameter; default = 2.403
#' 
#' probability of clinical infection/clinical immunity:
#'
#' * phi0 - maximum probability due to no immunity; default = 0.792
#' * phi1 - maximum reduction due to immunity; default = 0.00074
#' * ic0 - scale parameter; default = 18.02366
#' * kc - shape parameter; default = 2.36949
#'
#' probability of severe infection/severe disease immunity (P.f only):
#'
#' * theta0 - maximum probability due to no immunity; default = 0.0749886
#' * theta1 - maximum reduction due to immunity; default = 0.0001191
#' * iv0 - scale parameter; default = 1.09629
#' * kv - shape parameter; default = 2.00048
#' * fv0 - age dependent modifier; default = 0.141195
#' * av - age dependent modifier; default = 2493.41
#' * gammav - age dependent modifier; default = 2.91282
#'
#' infectivity towards mosquitoes:
#'
#' * cd - infectivity of clinically diseased humans towards mosquitoes; default = 0.068
#' * ca - infectivity of asymptomatically diseased humans towards mosquitoes (P.v only); default = 0.1
#' * gamma1 - parameter for infectivity of asymptomatic humans (P.f only); default = 1.82425
#' * cu - infectivity of sub-patent infection; default = 0.0062
#' * ct - infectivity of treated infection; default = 0.021896
#'
#' parasite incubation periods:
#'
#' * de - Duration of the human latent period of infection; default = 12
#' * delay_gam - Lag from parasites to infectious gametocytes; default = 12.5
#' * dem - Extrinsic incubation period in mosquito population model; default = 10
#' 
#' hypnozoite parameters (P.v only):
#' 
#' * f - relapse rate; default = 0.024
#' * gammal - clearance rate; default = 0.0026
#' 
#' parasite parameter
#' 
#' * parasite - parasite species (falciparum or vivax); default = "falciparum"
#' 
#' human population parameters:
#' * human_population - the initial number of humans to model; default = 100
#' * human_population_timesteps - the timesteps at which the population should
#' change; default = 0
#' * custom_demography - population demography given; default = FALSE,
#' * average_age - the average age of humans (in timesteps), this is only used
#' if custom_demography is FALSE; default = 7663
#' 
#' initial state proportions:
#'
#' * s_proportion - the proportion of `human_population` that begin as susceptible; default = 0.420433246
#' * d_proportion - the proportion of `human_population` that begin with
#' clinical disease; default = 0.007215064
#' * a_proportion - the proportion of `human_population` that begin as
#' asymptomatic; default = 0.439323667
#' * u_proportion - the proportion of `human_population` that begin as
#' subpatents; default = 0.133028023
#' * t_proportion - the proportion of `human_population` that begin treated; default = 0
#'
#' rendering:
#' All values are in timesteps and all ranges are inclusive
#' please set rendered age groups using the convenience function
#' `set_demography`
#'
#' * age_group_rendering_min_ages - the minimum ages for population size outputs; default = numeric(0)
#' * age_group_rendering_max_ages - the corresponding max ages; default = numeric(0)
#' * prevalence_rendering_min_ages - the minimum ages for clinical prevalence
#' outputs; default = 730
#' * prevalence_rendering_max_ages - the corresponding max ages; default = 3650
#' * incidence_rendering_min_ages - the minimum ages for incidence
#' outputs (includes asymptomatic microscopy +); default = turned off
#' * incidence_rendering_max_ages - the corresponding max ages; default = turned off
#' * clinical_incidence_rendering_min_ages - the minimum ages for clinical incidence outputs (symptomatic); default = 0
#' * clinical_incidence_rendering_max_ages - the corresponding max ages; default = 1825
#' * severe_incidence_rendering_min_ages - the minimum ages for severe incidence
#' outputs; default = turned off
#' * severe_incidence_rendering_max_ages - the corresponding max ages; default = turned off
#' * severe_prevalence_rendering_min_ages - the minimum ages for severe prevalance outputs; default = numeric(0),
#' * severe_prevalence_rendering_max_ages - the corresponding max ages; default = numeric(0),
#'
#' mosquito life stage transitions:
#'
#' * del - the delay for mosquitoes to move from state E to L; default = 6.64
#' * dl - the delay for mosquitoes to move from state L to P; default = 3.72
#' * dpl - the delay mosquitoes to move from state P to Sm; default = 0.643
#' * mup - the rate at which pupal mosquitoes die; default = 0.249
#' * mum - the rate at which developed mosquitoes die; default = 0.132
#'
#' unique biting rate:
#'
#' * a0 - age dependent biting parameter; default = 2920
#' * rho - age dependent biting parameter; default = 0.85
#' * sigma_squared - heterogeneity parameter; default = 1.67
#' * n_heterogeneity_groups - number discretised groups for heterogeneity, used
#' for sampling mothers; default = 5
#'
#' seasonality and carrying capacity parameters:
#'
#' * model_seasonality - boolean switch TRUE iff the simulation models seasonal rainfall; default = FALSE
#' * g0 - rainfall fourier parameter; default = 2
#' * g - rainfall fourier parameter; default = 0.3, 0.6, 0.9
#' * h - rainfall fourier parameters; default = 0.1, 0.4, 0.7
#' * gamma - effect of density dependence on late instars relative to early
#' instars; default = 13.25
#' * rainfall_floor - the minimum rainfall value (must be above 0); default 0.001
#'
#' flexible carrying capacity parameters:
#' please set carrying capacity using the convenience function
#' `set_carrying_capacity`
#' 
#' * carrying_capacity; default = FALSE
#' * carrying_capacity_timesteps; default = NULL
#' * carrying_capacity_values; default = NULL
#'  
#' larval mortality rates:
#'
#' * me - early stage larval mortality rate; default = 0.0338
#' * ml - late stage larval mortality rate; default = 0.0348
#'
#' vector biology:
#' species specific values are vectors
#' please set species parameters using the convenience function
#' `set_species`
#'
#' * beta - the average number of eggs laid per female mosquito per day; default = 21.2
#' * total_M - the initial number of adult mosquitos in the simulation; default = 1000
#' * init_foim - the FOIM used to calculate the equilibrium state for mosquitoes; default = 0
#' * species - names of the species in the simulation; default = "All"
#' * species_proportions - the relative proportions of each species; default = 1
#' * blood_meal_rates - the blood meal rates for each species; default = 0.3333333333
#' * Q0 - proportion of blood meals taken on humans; default = 0.92
#' * foraging_time - time spent taking blood meals; default = 0.69
#'
#' feeding cycle:
#' please set vector control strategies using `set_betnets` and `set_spraying`
#'
#' * bednets - boolean for if bednets are enabled; default = FALSE
#' * phi_bednets - proportion of bites taken in bed; default = 0.85
#' * k0 - proportion of females bloodfed with no net; default = 0.699
#' * spraying - boolean for if indoor spraying is enabled; default = FALSE
#' * phi_indoors - proportion of bites taken indoors; default = 0.90
#'
#' treatment parameters:
#' please set drug and treatment parameters with the convenience functions
#' `set_drugs` and `set_clinical_treatment`
#'
#' * drug_efficacy - a vector of efficacies for available drugs; default = turned off
#' * drug_rel_c - a vector of relative onward infectiousness values for drugs; default = turned off
#' * drug_prophylaxis_shape - a vector of shape parameters for weibull curves to
#' model prophylaxis for each drug; default = turned off
#' * drug_prophylaxis_scale - a vector of scale parameters for weibull curves to
#' model prophylaxis for each drug; default = turned off
#' * clinical_treatment_drugs - a list of drugs that are available for
#' clinically diseased (these values refer to the index in drug_* parameters); default = NULL
#' * clinical_treatment_timesteps - a list of vectors giving timesteps at which the
#' clinical treatment coverage changes (these values refer to the index in drug_* parameters); default = NULL
#' * clinical_treatment_coverage - a list of vectors giving coverage values for each drug; default = NULL
#'
#' PEV parameters: 
#' please set vaccine strategies with the convenience functions
#' `set_pev_epi` and `set_mass_pev`
#'
#' * pev - PEV strategy implemented; default = FALSE
#' * pev_doses - the dosing schedule before the vaccine takes effect; default =
#' c(0, 1.5 * 30, 3 * 30)
#' default = 365
#'
#' MDA, SMC and PMC parameters:
#' please set these parameters with the convenience functions
#' `set_mda`, `set_smc` and `set_pmc`,
#' with `peak_season_offset `
#'
#' TBV parameters:
#' please set TBV parameters with the convenience functions
#' `set_tbv`
#'
#' * tbv; default = FALSE
#' * tbv_mt - effect on treated infectiousness; default = 35
#' * tbv_md - effect on diseased infectiousness; default = 46.7
#' * tbv_ma - effect on asymptomatic infectiousness; default = 3.6
#' * tbv_mu - effect on subpatent infectiousness; default = 0.8
#' * tbv_k  - scale parameter for effect on infectiousness; default = 0.9
#' * tbv_tau - peak antibody parameter; default = 22
#' * tbv_rho - antibody component parameter; default = 0.7
#' * tbv_ds - antibody short-term delay parameter; default = 45
#' * tbv_dl - antibody long-term delay parameter; default = 591
#' * tbv_tra_mu - transmission reduction parameter; default = 12.63
#' * tbv_gamma1 - transmission reduction parameter; default = 2.5
#' * tbv_gamma2 - transmission reduction parameter; default = 0.06
#' * tbv_timesteps; default = NULL
#' * tbv_coverages; default = NULL
#' * tbv_ages; default = NULL
#'
#'
#' miscellaneous:
#' please set parasite species using the convenience function
#' `set_parasite`
#' 
#' * mosquito_limit - the maximum number of mosquitoes to allow for in the
#' simulation; default = 1.00E+05
#' * individual_mosquitoes - boolean whether adult mosquitoes are modeled
#' individually or compartmentally; default = TRUE
#' * enable_heterogeneity - boolean whether to include heterogeneity in biting
#' rates; default = TRUE
#' * r_tol - the relative tolerance for the ode solver; default = 1e-4
#' * a_tol - the absolute tolerance for the ode solver; default = 1e-4
#' * ode_max_steps - the max number of steps for the solver; default = 1e6
#' *  progress_bar - display progress bar; default = FALSE
#'
#' @export
get_parameters <- function(overrides = list(), parasite = "falciparum") {
  
  parameters <- c(
    ## Parasite-specific parameters set in parasite_parameters.csv
    # human fixed state transitions
    # duration of sub-patent infection
    # initial immunity values
    # immunity boost grace periods
    # maternal immunity parameters
    # immunity decay rates
    # probability of pre-erythrocytic infection/blood immunity
    # probability of asymptomatic detection/detectbility immunity 
    # probability of asymptomatic infection/antiparasite immunity
    # probability of clinical infection/clinical immunity
    # probability of severe infection/severe disease immunity
    # infectivity towards mosquitos
    # parasite incubation periods
    # hypnozoite parameters
    parasite_parameters[[parasite]],
    
    ## Add remaining parameters
    list(
      # parasite parameter
      parasite = parasite,
      
      # human population parameters
      human_population = 100,
      human_population_timesteps = 0,
      custom_demography = FALSE,
      average_age = 7663,
      # initial state proportions
      s_proportion = 0.420433246,
      d_proportion = 0.007215064,
      a_proportion = 0.439323667,
      u_proportion = 0.133028023,
      t_proportion = 0,
      # rendering
      age_group_rendering_min_ages = numeric(0),
      age_group_rendering_max_ages = numeric(0),
      prevalence_rendering_min_ages = 2 * 365,
      prevalence_rendering_max_ages = 10 * 365,
      incidence_rendering_min_ages = numeric(0),
      incidence_rendering_max_ages = numeric(0),
      clinical_incidence_rendering_min_ages = numeric(0),
      clinical_incidence_rendering_max_ages = 5 * 365,
      severe_incidence_rendering_min_ages = numeric(0),
      severe_incidence_rendering_max_ages = numeric(0),
      severe_prevalence_rendering_min_ages = numeric(0),
      severe_prevalence_rendering_max_ages = numeric(0),
      
      # mosquito life stage transitions
      del =  6.64,
      dl  =  3.72,
      dpl =  0.643,
      mup =  0.249,
      mum =  0.132,
      # unique biting rate
      a0  =  2920,
      rho =  0.85,
      sigma_squared  = 1.67,
      n_heterogeneity_groups  =  5,
      # seasonality and carrying capacity parameters
      model_seasonality = FALSE,
      g0    = 2,
      g     = c(.3, .6, .9),
      h     = c(.1, .4, .7),
      gamma = 13.25,
      rainfall_floor = 0.001,
      # flexible carrying capacity parameters
      carrying_capacity = FALSE,
      carrying_capacity_timesteps = NULL,
      carrying_capacity_values = NULL,
      # larval mortality rates
      me    = .0338,
      ml    = .0348,
      # vector biology
      beta     = 21.2,
      total_M  = 1000,
      init_foim= 0,
      # order of species: An gambiae s.s, An arabiensis, An funestus
      species             = 'All',
      species_proportions = 1,
      blood_meal_rates    = 1/3,
      Q0                  = .92,
      foraging_time       = .69,
      # bed nets
      bednets = FALSE,
      phi_bednets = .85,
      k0 = .699,
      # indoor spraying
      spraying = FALSE,
      phi_indoors = .90,
      # treatment
      drug_efficacy          = numeric(0),
      drug_rel_c             = numeric(0),
      drug_prophylaxis_shape = numeric(0),
      drug_prophylaxis_scale = numeric(0),
      clinical_treatment_drugs     = list(),
      clinical_treatment_timesteps = list(),
      clinical_treatment_coverages = list(),
      # rts,s
      pev = FALSE,
      pev_doses = c(0, 1.5 * 30, 3 * 30),
      # MDA
      mda = FALSE,
      mda_drug = 0,
      mda_timesteps = NULL,
      mda_coverages = NULL,
      mda_min_ages = -1,
      mda_max_ages = -1,
      smc = FALSE,
      smc_drug = 0,
      smc_timesteps = NULL,
      smc_coverages = NULL,
      smc_min_ages = -1,
      smc_max_ages = -1,
      # PMC
      pmc = FALSE,
      pmc_drug = 0,
      pmc_timesteps = NULL,
      pmc_coverages = NULL,
      pcs_ages = -1,
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
      tbv_timesteps = NULL,
      tbv_coverages = NULL,
      tbv_ages = NULL,
      
      # misc
      mosquito_limit   = 100 * 1000,
      individual_mosquitoes = FALSE,
      enable_heterogeneity = TRUE,
      r_tol = 1e-4,
      a_tol = 1e-4,
      ode_max_steps = 1e6,
      progress_bar = FALSE
      ))
  
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
  parameterise_total_M(parameters, equilibrium_total_M(parameters, EIR))
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
  if (!parameters$individual_mosquitoes) {
    return(parameters)
  }
  max_total_M <- 0
  for (i in seq_along(parameters$species)) {
    species_M <- total_M * parameters$species_proportions[[i]]
    K0 <- calculate_carrying_capacity(parameters, species_M, i)
    R_bar <- calculate_R_bar(parameters)
    max_K <- max(vnapply(seq(365), function(t) {
      carrying_capacity(
        t,
        parameters$model_seasonality,
        parameters$g0,
        parameters$g,
        parameters$h,
        K0,
        R_bar,
        parameters$rainfall_floor
      )
    }))
    omega <- calculate_omega(parameters, i)
    mum <- weighted.mean(parameters$mum, parameters$species_proportions)
    max_total_M <- max_total_M + max_K * (
      1 / (
        2 * parameters$dl * mum * (
          1 + parameters$dpl * parameters$mup
        )
      )
    ) * (
      1 / (
        parameters$gamma * (omega + 1)
      )
    ) * (
      omega / (parameters$ml * parameters$del) - (
        1 / (parameters$ml * parameters$dl)
      ) - 1
    )
  }
  parameters$mosquito_limit <- ceiling(max_total_M * 5) #Allow for random fluctuations
  parameters
}

#' Use parameter draw from the join posterior
#' 
#' Overrides default (median) model parameters with a single draw from the fitted
#' joint posterior. Must be called prior to set_equilibrium.
#'
#' @param parameters the model parameters
#' @param draw the draw to use. Must be an integer between 1 and 1000
#'
#' @export
set_parameter_draw <- function(parameters, draw){
  if(draw > 1000 || draw < 1){
    stop("draw must be an integer between 1 and 1000")
  }
  parameter_draw <- parameter_draws[[draw]]
  for (name in names(parameter_draw)) {
    parameters[[name]] <- parameter_draw[[name]]
  }
  return(parameters)
}


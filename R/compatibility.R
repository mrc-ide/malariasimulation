EQUILIBRIUM_AGES <- 0:999 / 10

inverse_param <- function(name, new_name) {
  function(params) { list(new_name, 1 / params[[name]]) }
}

product_param <- function(new_name, name_1, name_2) {
  function(params) { list(new_name, params[[name_1]] * params[[name_2]]) }
}

mean_param <- function(new_name, name, weights) {
  function(params) {
    list(new_name, weighted.mean(params[[name]], params[[weights]]))
  }
}

translations = list(
  eta = inverse_param('eta', 'average_age'),
  rho = 'rho',
  a0  = 'a0',
  rA  = inverse_param('rA', 'da'),
  rD  = inverse_param('rD', 'dd'),
  rU  = inverse_param('rU', 'du'),
  rT  = inverse_param('rT', 'dt'),
  dE  = 'de',
  cD  = 'cd',
  cU  = 'cu',
  cT  = 'ct',
  d1  = 'd1',
  dd  = 'rid',
  ID0 = 'id0',
  kd  = 'kd',
  ud  = 'ud',
  ad0 = 'ad',
  gd  = 'gammad',
  b0  = 'b0',
  b1  = 'b1',
  d1  = 'd1',
  db  = 'rb',
  IB0 = 'ib0',
  kb  = 'kb',
  ub  = 'ub',
  phi0= 'phi0',
  phi1= 'phi1',
  dc  = 'rc',
  IC0 = 'ic0',
  kc  = 'kc',
  uc  = 'uc',
  dm  = 'rm',
  mu  = 'mum',
  s2  = 'sigma_squared',
  fd0 = 'fd0',
  g_inf = 'gamma1',
  PM = 'pcm',
  tau = 'dem',
  Q0 = 'Q0',
  f = 'blood_meal_rates',
  tl = 'delay_gam',
  uv = 'uv',
  dv = 'rva',
  PVM = 'pvm',
  dvm = 'rvm',
  fv0 = 'fv0',
  av = 'av',
  gammav = 'gammav',
  theta0 = 'theta0',
  theta1 = 'theta1',
  IV0 = 'iv0',
  kv = 'kv'
)

back_translations = list(
  average_age = inverse_param('average_age', 'eta'),
  rho = 'rho',
  a0  = 'a0',
  da = inverse_param('da', 'rA'),
  dd = inverse_param('dd', 'rD'),
  du  = inverse_param('du', 'rU'),
  dt  = inverse_param('dt', 'rT'),
  de  = 'dE',
  cd  = 'cD',
  cu  = 'cU',
  ct  = 'cT',
  d1  = 'd1',
  rid = 'dd',
  id0 = 'ID0',
  kd  = 'kd',
  ud  = 'ud',
  ad  = 'ad0',
  gammad  = 'gd',
  b0  = 'b0',
  b1  = 'b1',
  d1  = 'd1',
  rb  = 'db',
  ib0 = 'IB0',
  kb  = 'kb',
  ub  = 'ub',
  phi0= 'phi0',
  phi1= 'phi1',
  rc  = 'dc',
  ic0 = 'IC0',
  kc  = 'kc',
  uc  = 'uc',
  rm  = 'dm',
  mum  = mean_param('mu', 'mum', 'species_proportions'),
  sigma_squared = 's2',
  fd0 = 'fd0',
  gamma1 = 'g_inf',
  pcm = 'PM',
  dem = 'tau',
  Q0 = mean_param('Q0', 'Q0', 'species_proportions'),
  blood_meal_rates = mean_param('f', 'blood_meal_rates', 'species_proportions'),
  delay_gam = 'tl',
  uv = 'uv',
  rva = 'dv',
  pvm = 'PVM',
  rvm = 'dvm',
  fv0 = 'fv0',
  av = 'av',
  gammav = 'gammav',
  theta0 = 'theta0',
  theta1 = 'theta1',
  IV0 = 'iv0',
  kv = 'kv'
)

vivax_translations = list(
  
  mean_age  = 'average_age',
  rho_age   = 'rho',
  age_0     = 'a0',
  N_het     = 'n_heterogeneity_groups',
  
  bb  = 'b',      ## mosquito -> human transmission probability
  
  c_PCR = 'cu',   ## human -> mosquito transmission probability (PCR)
  c_LM  = 'ca',   ## human -> mosquito transmission probability (LM-detectable)
  c_D   = 'cd',   ## human -> mosquito transmission probability (disease state)
  c_T   = 'ct',   ## human -> mosquito transmission probability (treatment)
  
  d_E   = 'de',   ## duration of liver-stage latency
  r_D   = inverse_param('dd', 'r_D'),   ## duraton of disease = 1/rate
  r_T   = inverse_param('dt', 'r_T'),   ## duraton of prophylaxis = 1/rate
  
  r_par = inverse_param('ra', 'r_par'),    ## rate of decay of anti-parasite immunity
  r_clin= inverse_param('rc', 'r_clin'),    ## rate of decay of clinical immunity
  
  mu_M  = mean_param('mum', 'mum', 'species_proportions'),  ## mosquito death rate = 1/(mosquito life expectancy)
  Q0    = mean_param('Q0', 'Q0', 'species_proportions'),
  blood_meal_rates = mean_param('blood_meal_rates', 'blood_meal_rates', 'species_proportions'),
  tau_M = 'dem', 	## duration of sporogony
  
  ff      = 'f',      ## relapse rate
  gamma_L = 'gammal', ## duration of liver-stage carriage
  K_max   = 'kmax', ## maximum hypnozoite batches
  
  u_par        = 'ua',        ## refractory period for anti-parasite immune boosting
  phi_LM_max   = 'philm_max', ## probability of LM_detectable infection with no immunity
  phi_LM_min   = 'philm_min', ## probability of LM_detectable infection with maximum immunity
  A_LM_50pc    = 'alm50',     ## blood-stage immunity scale parameter
  K_LM         = 'klm',       ## blood-stage immunity shape parameter
  u_clin       = 'uc',        ## refractory period for clinical immune boosting
  phi_D_max    = 'phi0',      ## probability of clinical episode with no immunity
  phi_D_min    = product_param("phi_D_min", "phi0", "phi1"),    ## probability of clinical episode with maximum immunity
  A_D_50pc     = 'ic0',       ## clinical immunity scale parameter
  K_D          = 'kc',        ## clinical immunity shape parameter
  A_d_PCR_50pc = 'apcr50',    ## scale parameter for effect of anti-parasite immunity on PCR-detectable infection
  K_d_PCR      = 'kpcr',      ## shape parameter for effect of anti-parasite immunity on PCR-detectable infection
  d_PCR_max    = 'dpcr_max',  ## maximum duration on PCR-detectable infection
  d_PCR_min    = 'dpcr_min',  ## maximum duration of PCR-detectable infection
  d_LM         = 'da',        ## duration of LM-detectable infection
  P_MI         = 'pcm',       ## Proportion of immunity acquired maternally
  d_MI         = 'rm'         ## Rate of waning of maternal immunity
  
)

#' @description translate parameter keys from the malariaEquilibrium format
#' to ones compatible with this IBM 
#' @param params with keys in the malariaEquilibrium format
#' @noRd
translate_equilibrium <- function(params) {
  translated <- list()
  for (name in names(params)) {
    if(!name %in% names(translations)) {
      stop(paste('Unknown parameter', name))
    }
    translation <- translations[[name]]
    if (is.character(translation)) {
      translated[[translation]] <- params[[name]]
    }
    if (is.function(translation)) {
      t <- translation(params)
      translated[[t[[1]]]] <- t[[2]]
    }
  }
  translated
}

#' @description translate model parameters to the malariaEquilibrium format
#' @param params model params
#' @noRd
translate_parameters <- function(params) {
  translated <- malariaEquilibrium::load_parameter_set()
  for (name in names(params)) {
    if(name %in% names(back_translations)) {
      translation <- back_translations[[name]]
      if (is.character(translation)) {
        translated[[translation]] <- params[[name]]
      }
      if (is.function(translation)) {
        t <- translation(params)
        translated[[t[[1]]]] <- t[[2]]
      }
    }
  }
  translated
}

#' @description translate parameter keys from the malariaVivaxEquilibrium format
#' to ones compatible with this IBM
#' @param params with keys in the malariaVivaxEquilibrium format
#' @noRd
translate_vivax_parameters <- function(params) {
  translated <- params
  for (i in 1:length(vivax_translations)) {
    if (is.character(vivax_translations[[i]])) {
      translated[[names(vivax_translations)[i]]] <- params[[vivax_translations[[i]]]]
    }
    if (is.function(vivax_translations[[i]])) {
      translated[[names(vivax_translations)[i]]] <- vivax_translations[[i]](params)[[2]]
    }
  }
  translated
}

#' @title remove parameter keys from the malariaEquilibrium format that are not used
#' in this IBM 
#' @param params with keys in the malariaEquilibrium format
#' @noRd
remove_unused_equilibrium <- function(params) {
  remove_keys(
    params,
    c(
      'rP', # Prophylaxis state is no longer used, see `drug_parameters.R`
      'aA', # used for pcr calculations
      'aU', # used for pcr calculations
      'cd_w', # unused!
      'cd_p'  # unused!
    )
  )
}

#' @title Set equilibrium
#' @description This will update the IBM parameters to match the
#' equilibrium parameters and set up the initial human and mosquito population
#' to achieve init_EIR
#' @param parameters model parameters to update
#' @param init_EIR the desired initial EIR (infectious bites per person per day over the entire human
#' population)
#' @param eq_params parameters from the malariaEquilibrium package, if null.
#' The default malariaEquilibrium parameters will be used
#' @export
set_equilibrium <- function(parameters, init_EIR, eq_params = NULL) {
  if(parameters$parasite == "falciparum"){
    if (is.null(eq_params)) {
      eq_params <- translate_parameters(parameters)
    } else {
      parameters <- c(
        translate_equilibrium(remove_unused_equilibrium(eq_params)),
        parameters
      )
    }
    eq <- malariaEquilibrium::human_equilibrium(
      EIR = init_EIR,
      ft = sum(get_treatment_coverages(parameters, 1)),
      p = eq_params,
      age = EQUILIBRIUM_AGES,
      h = malariaEquilibrium::gq_normal(parameters$n_heterogeneity_groups)
    )
    parameters <- c(
      list(
        init_foim = eq$FOIM,
        init_EIR = init_EIR,
        eq_params = eq_params
      ),
      parameters
    )
  } else if (parameters$parasite == "vivax"){
    
    if (!(is.null(eq_params))) {
      stop("Importing MalariaEquilibriumVivax parameters is not supported")
    }
    
    eq <- malariaEquilibriumVivax::vivax_equilibrium(
      EIR = init_EIR,
      ft = sum(get_treatment_coverages(parameters, 1)),
      p = translate_vivax_parameters(parameters),
      age = EQUILIBRIUM_AGES
    )
    
    parameters <- c(
      list(
        init_foim = eq$FOIM,
        init_EIR = init_EIR
      ),
      parameters
    )
  }
  parameterise_mosquito_equilibrium(parameters, init_EIR)
}

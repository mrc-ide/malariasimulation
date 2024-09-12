EQUILIBRIUM_AGES <- 0:999 / 10

inverse_param <- function(name, new_name) {
  function(params) { list(new_name, 1 / params[[name]]) }
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
  blood_meal_rates = mean_param('f', 'blood_meal_rates', 'semiochemical_effect', 'species_proportions'),
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
#' to acheive init_EIR
#' @param parameters model parameters to update
#' @param init_EIR the desired initial EIR (infectious bites per person per day over the entire human
#' population)
#' @param eq_params parameters from the malariaEquilibrium package, if null.
#' The default malariaEquilibrium parameters will be used
#' @export
set_equilibrium <- function(parameters, init_EIR, eq_params = NULL) {
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
  parameterise_mosquito_equilibrium(parameters, init_EIR)
}

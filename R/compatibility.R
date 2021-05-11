EQUILIBRIUM_AGES <- 0:999 / 10

inverse_param <- function(name, new_name) {
  function(params) { list(new_name, 1 / params[[name]]) }
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
  tl = 'delay_gam'
)

#' @title translate parameter keys from the malariaEquilibrium format to ones compatible
#' with this IBM 
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
#' @param init_EIR the desired initial EIR
#' @param eq_params parameters from the malariaEquilibrium package, if null.
#' The default malariaEquilibrium parameters will be used
#' @export
set_equilibrium <- function(parameters, init_EIR, eq_params = NULL) {
  if (is.null(eq_params)) {
    eq_params <- malariaEquilibrium::load_parameter_set()
  }
  eq <- malariaEquilibrium::human_equilibrium(
    EIR = init_EIR,
    ft = sum(get_treatment_coverages(parameters, 1)),
    p = eq_params,
    age = EQUILIBRIUM_AGES,
    h = malariaEquilibrium::gq_normal(parameters$n_heterogeneity_groups)
  )
  parameters <- c(
    translate_equilibrium(remove_unused_equilibrium(eq_params)),
    list(
      init_foim = eq$FOIM,
      init_EIR = init_EIR,
      eq_params = eq_params
    ),
    parameters
  )
  parameterise_mosquito_equilibrium(parameters, init_EIR)
}

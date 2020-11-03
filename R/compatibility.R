
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
  f = 'blood_meal_rates'
)

#' @title translate parameter keys from jamie's format to ones compatible
#' with this IBM 
#' @param params with keys in the jamie's format
#' @export
translate_jamie <- function(params) {
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

#' @title remove parameter keys from jamie's format that are not used
#' in this IBM 
#' @param params with keys in the jamie's format
#' @export
remove_unused_jamie <- function(params) {
  remove_keys(
    params,
    c(
      'rP', # Prophylaxis state is no longer used, see `drug_parameters.R`
      'tl', # unused!
      'aA', # used for microscopy and pcr calculations
      'aU', # used for microscopy and pcr calculations
      'cd_w', # unused!
      'cd_p'  # unused!
    )
  )
}

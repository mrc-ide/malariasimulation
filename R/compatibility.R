
inverse_param <- function(name, new_name) {
  function(params) { list(new_name, 1 / params[[name]]) }
}

translations = list(
  eta = 'human_death_rate',
  rho = 'rho',
  a0  = 'a0',
  rA  = inverse_param('rA', 'da'),
  rD  = inverse_param('rD', 'dd'),
  rU  = inverse_param('rU', 'du'),
  dE  = 'rel',
  cD  = 'cd',
  cU  = 'cu',
  d1  = 'd1',
  dd  = inverse_param('dd', 'rid'),
  ID0 = 'id0',
  kd  = 'kd',
  ud  = 'ud',
  ad0 = 'ad',
  gd  = 'gammad',
  b0  = 'b0',
  d1  = 'd1',
  db  = inverse_param('db', 'rb'),
  IB0 = 'id0',
  kb  = 'kb',
  ub  = 'ub',
  phi0= 'phi0',
  phi1= 'phi1',
  dc  = inverse_param('dc', 'rc'),
  IC0 = 'ic0',
  kc  = 'kc',
  uc  = 'uc',
  dm  = inverse_param('dm', 'rm'),
  mu  = 'mum'
)

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

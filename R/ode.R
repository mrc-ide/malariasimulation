parameterise_ode <- function(parameters, foim = 0.) {
  create_mosquito_model(
    initial_mosquito_counts(parameters, foim),
    parameters$beta,
    parameters$del,
    parameters$me,
    parameters$K0,
    parameters$gamma,
    parameters$dl,
    parameters$ml,
    parameters$dpl,
    parameters$mup,
    foim,
    parameters$mum,
    parameters$dem
  )
}

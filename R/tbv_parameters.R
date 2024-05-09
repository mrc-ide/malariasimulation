#' @title Parameterise an TBV strategy
#' @param parameters a list of parameters to modify
#' @param timesteps a vector of timesteps for each round of vaccinations
#' @param coverages the coverage for each round of vaccinations
#' @param ages a vector of ages of the target population (in years)
#' @export
set_tbv <- function(
  parameters,
  timesteps,
  coverages,
  ages
  ) {
  stopifnot(all(coverages >= 0) && all(coverages <= 1))
  if(length(coverages) != length(timesteps)){
    stop("coverages and timesteps do no align")
  }

  parameters$tbv <- TRUE
  parameters$tbv_timesteps <- timesteps
  parameters$tbv_coverages <- coverages
  parameters$tbv_ages <- ages
  parameters
}

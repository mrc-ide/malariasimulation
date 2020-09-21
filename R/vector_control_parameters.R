#' @title Parameterise a bed net strategy
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to distribute bednets
#' @param coverages the proportion of the target population who receive bednets
#' @param retention the average number of timesteps a net is kept for
#' @export
set_bednets <- function(
  parameters,
  timesteps,
  coverages,
  retention
  ) {
  if (length(timesteps) != length(coverages)) {
    stop('timesteps and coverages must align')
  }
  parameters$bednets <- TRUE
  parameters$bednet_timesteps <- timesteps
  parameters$bednet_coverages <- coverages
  parameters$bednet_retention <- retention
  parameters
}

#' @title Parameterise an indoor spraying strategy
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to spray
#' @param coverages the proportion of the target population who get indoor
#' spraying
#' @export
set_spraying <- function(
  parameters,
  timesteps,
  coverages
  ) {
  if (length(timesteps) != length(coverages)) {
    stop('timesteps and coverages must align')
  }
  parameters$spraying <- TRUE
  parameters$spraying_timesteps <- timesteps
  parameters$spraying_coverages <- coverages
  parameters
}

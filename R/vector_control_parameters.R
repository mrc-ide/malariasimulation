#' @title Parameterise a bed net strategy
#'
#' @description The model will distribute bednets at `timesteps` to a random
#' sample of the entire human population. The sample size will be a proportion
#' of the human population taken from the corresponding `coverages`.
#' The sample _can_ contain humans who already have bednets.
#'
#' All of the sample "use" their bednets on the timestep after they are
#' distributed. Incomplete usage is not part of this model.
#'
#' If a human in the sample already has a bednet, their bednet will be replaced
#' by a new one.
#'
#' Bednets will be randomly removed each timestep with a rate of `1 -
#' exp(-1/retention)`
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to distribute bednets
#' @param coverages the proportion of the human population who receive bednets
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
#'
#' @description The model will apply indoor spraying at `timesteps` to a random
#' sample of the entire human population. The sample size will be a proportion
#' of the human population taken from the corresponding `coverages`.
#' The sample _can_ contain humans who have already benefitted from spraying.
#'
#' If a human in the sample lives in a sprayed house, the efficacy of the
#' spraying will be returned to the maximum.
#'
#' @param parameters a list of parameters to modify
#' @param timesteps the timesteps at which to spray
#' @param coverages the proportion of the population who get indoor
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

#' @title Parameterise an RTS,S strategy
#' @param parameters a list of parameters to modify
#' @param start the timestep to start
#' @param end the last timestep for the intervention
#' @param frequency the number of timsteps between rounds
#' @param ages for each round (in years)
#' @param coverage the proportion of the target population who will recieve the
#' vaccine
#' @export
set_rtss <- function(
  parameters,
  start,
  end,
  frequency,
  ages,
  coverage
  ) {
  parameters$rtss <- TRUE
  parameters$rtss_start <- start
  parameters$rtss_end <- end
  parameters$rtss_frequency <- frequency
  parameters$rtss_ages <- ages
  parameters$rtss_coverage <- coverage
  parameters
}

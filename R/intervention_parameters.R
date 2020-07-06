parameter_push_back <- function(parameters, name, value) {
  parameters[[name]] <- c(parameters[[name]], value)
  parameters
}

#' @title Parameterise a drug for use in human interventions
#' @param parameters a list of parameters to modify
#' @param efficacy of the drug to add
#' @param efficacy of the drug to add
#' @param infectiousness of the individual during treatment
#' @export
add_drug <- function(parameters, efficacy) {
  parameter_push_back(parameters, 'drug_efficacies', efficacy)
}

#' @title Parameterise a Mass Drug Administration
#' @param parameters a list of parameters to modify
#' @param drug the index of the drug to administer
#' @param start the timestep to start
#' @param end the last timestep for the intervention
#' @param frequency the number of timsteps between doses
#' @param min_age the minimum age of the target population exclusive (in timesteps)
#' @param max_age the maximum age of the target population exclusive (in timesteps)
#' @param coverage the proportion of the target population who will recieve the
#' drug
#' @export
add_mda <- function(
  parameters,
  drug,
  start,
  end,
  frequency,
  min_age,
  max_age,
  coverage
  ) {
  parameters <- parameter_push_back(parameters, 'mda_drug', drug)
  parameters <- parameter_push_back(parameters, 'mda_start', start)
  parameters <- parameter_push_back(parameters, 'mda_end', end)
  parameters <- parameter_push_back(parameters, 'mda_frequency', frequency)
  parameters <- parameter_push_back(parameters, 'mda_min_age', min_age)
  parameters <- parameter_push_back(parameters, 'mda_max_age', max_age)
  parameter_push_back(parameters, 'mda_coverage', coverage)
}

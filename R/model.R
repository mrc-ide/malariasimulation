#' @title Run the simulation
#' @description
#' Run the simulation for some time given some parameters. This currently
#' returns a dataframe with the number of individuals in each state at each
#' timestep.
#'
#' The resulting dataframe contains the following columns:
#'
#'  * timestep: the timestep for the row
#'  * infectivity: the infectivity from humans towards mosquitoes
#'  * lambda: the effective biting rate on humans (per timestep) (per
#'species). This defaults to lambda_All
#'  * normal_lambda: the effective biting rate on adult humans (per
#' timestep) (per species). This defaults to normal_lambda_All
#'  * FOIM: the force of infection towards mosquitoes (per species)
#'  * mu: the death rate of adult mosquitoes (per species)
#'  * EIR: the Entomological Inoculation Rate (per timestep, over the whole
#'  population)
#'  * n_bitten: number of humans bitten by an infectious mosquito
#'  * n_treated: number of humans treated for clinical or severe malaria this timestep
#'  * n_infections: number of humans who get an asymptomatic, clinical or severe malaria this timestep
#'  * natural_deaths: number of humans who die from aging
#'  * S_count: number of humans who are Susceptible
#'  * A_count: number of humans who are Asymptomatic
#'  * D_count: number of humans who have the clinical malaria
#'  * U_count: number of subpatent infections in humans
#'  * Tr_count: number of infections being treated in humans
#'  * ica_mean: the mean acquired immunity to clinical infection over the population of humans
#'  * icm_mean: the mean maternal immunity to clinical infection over the population of humans
#'  * ib_mean: the mean blood immunity to all infection over the population of humans
#'  * id_mean: the mean immunity from detection through microscopy over the population of humans
#'  * n: number of humans between an inclusive age range at this timestep. This
#' defaults to n_730_3650. Other age ranges can be set with
#' prevalence_rendering_min_ages and prevalence_rendering_max_ages parameters.
#' * n_detect: number of humans with an infection detectable by microscopy between an inclusive age range at this timestep. This
#' defaults to n_detect_730_3650. Other age ranges can be set with
#' prevalence_rendering_min_ages and prevalence_rendering_max_ages parameters.
#'  * n_severe: number of humans with a severe infection detectable by microscopy 
#'  between an inclusive age range at this timestep. Age ranges can be set with
#' severe_prevalence_rendering_min_ages and severe_prevalence_rendering_max_ages parameters.
#'  * n_inc: number of new infections for humans between an inclusive age range at this timestep.
#' incidence columns can be set with
#' incidence_rendering_min_ages and incidence_rendering_max_ages parameters.
#'  * n_inc_clinical: number of new clinical infections for humans between an inclusive age range at this timestep. 
#' clinical incidence columns can be set with
#' clinical_incidence_rendering_min_ages and clinical_incidence_rendering_max_ages parameters.
#'  * n_inc_severe: number of new severe infections for humans between an inclusive age range at this timestep.
#' severe incidence columns can be set with
#' severe_incidence_rendering_min_ages and severe_incidence_rendering_max_ages parameters.
#'  * severe_deaths: number of deaths due to severe malaria. severe_enabled must be
#' set to TRUE
#'  * E_count: number of mosquitoes in the early larval stage
#'  * L_count: number of mosquitoes in the late larval stage
#'  * P_count: number of mosquitoes in the pupal stage
#'  * Sm_count: number of adult female mosquitoes who are Susceptible
#'  * Pm_count: number of adult female mosquitoes who are incubating
#'  * Im_count: number of adult female mosquitoes who are infectious
#'  * total_M: number of adult female mosquitoes. Variables with suffixes total_M_# refer to the number of adult 
#' female mosquitoes from the different indicies of the set_species() vector
#'  * rate_D_A: rate that humans transition from clinical disease to
#' asymptomatic
#'  * rate_A_U: rate that humans transition from asymptomatic to
#' subpatent
#'  * rate_U_S: rate that humans transition from subpatent to
#' susceptible
#'  * mosquito_deaths: number of adult female mosquitoes who die this timestep
#'
#' @param timesteps the number of timesteps to run the simulation for (in days)
#' @param parameters a named list of parameters to use
#' @param correlations correlation parameters
#' @return dataframe of results
#' @export
run_simulation <- function(
  timesteps,
  parameters = NULL,
  correlations = NULL
  ) {
  random_seed(ceiling(runif(1) * .Machine$integer.max))
  if (is.null(parameters)) {
    parameters <- get_parameters()
  }
  if (is.null(correlations)) {
    correlations <- get_correlation_parameters(parameters)
  }
  variables <- create_variables(parameters)
  events <- create_events(parameters)
  initialise_events(events, variables, parameters)
  renderer <- individual::Render$new(timesteps)
  attach_event_listeners(
    events,
    variables,
    parameters,
    correlations,
    renderer
  )
  vector_models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(vector_models, parameters)
  individual::simulation_loop(
    processes = create_processes(
      renderer,
      variables,
      events,
      parameters,
      vector_models,
      solvers,
      correlations
    ),
    variables = variables,
    events = unlist(events),
    timesteps = timesteps
  )
  renderer$to_dataframe()
}

#' @title Run the simulation with repetitions
#'
#' @param timesteps the number of timesteps to run the simulation for
#' @param repetitions n times to run the simulation
#' @param overrides a named list of parameters to use instead of defaults
#' @param parallel execute runs in parallel
#' @export
run_simulation_with_repetitions <- function(
  timesteps,
  repetitions,
  overrides = list(),
  parallel = FALSE
  ) {
  if (parallel) {
    fapply <- parallel::mclapply
  } else {
    fapply <- lapply
  }
  dfs <- fapply(
    seq(repetitions),
    function(repetition) {
      df <- run_simulation(timesteps, overrides)
      df$repetition <- repetition
      df
    }
  )
  do.call("rbind", dfs)
}

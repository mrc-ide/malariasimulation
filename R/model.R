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
#'  * FOIM: the force of infection towards mosquitoes (per species)
#'  * mu: the death rate of adult mosquitoes (per species)
#'  * EIR: the Entomological Inoculation Rate (per timestep, per species, over 
#'  the whole population)
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
#'  * n_detect: number of humans with an infection detectable by microscopy between an inclusive age range at this timestep. This
#' defaults to n_detect_730_3650. Other age ranges can be set with
#' prevalence_rendering_min_ages and prevalence_rendering_max_ages parameters.
#'  * p_detect: the sum of probabilities of detection by microscopy between an
#' inclusive age range at this timestep. This
#' defaults to p_detect_730_3650. Other age ranges can be set with
#' prevalence_rendering_min_ages and prevalence_rendering_max_ages parameters.
#'  * n_severe: number of humans with a severe infection between an inclusive
#' age range at this timestep. Age ranges can be set with
#' severe_prevalence_rendering_min_ages and severe_prevalence_rendering_max_ages parameters.
#'  * n_inc: number of new infections for humans between an inclusive age range at this timestep.
#' incidence columns can be set with
#' incidence_rendering_min_ages and incidence_rendering_max_ages parameters.
#'  * p_inc: sum of probabilities of infection for humans between an inclusive age range at this timestep.
#' incidence columns can be set with
#' incidence_rendering_min_ages and incidence_rendering_max_ages parameters.
#'  * n_inc_clinical: number of new clinical infections for humans between an inclusive age range at this timestep. 
#' clinical incidence columns can be set with
#' clinical_incidence_rendering_min_ages and clinical_incidence_rendering_max_ages parameters.
#'  * p_inc_clinical: sub of probabilities of clinical infection for humans between an inclusive age range at this timestep. 
#' clinical incidence columns can be set with
#' clinical_incidence_rendering_min_ages and clinical_incidence_rendering_max_ages parameters.
#'  * n_inc_severe: number of new severe infections for humans between an inclusive age range at this timestep.
#' severe incidence columns can be set with
#' severe_incidence_rendering_min_ages and severe_incidence_rendering_max_ages parameters.
#'  * p_inc_severe: the sum of probabilities of severe infection for humans between an inclusive age range at this timestep.
#' severe incidence columns can be set with
#' severe_incidence_rendering_min_ages and severe_incidence_rendering_max_ages parameters.
#'  * E_count: number of mosquitoes in the early larval stage (per species)
#'  * L_count: number of mosquitoes in the late larval stage (per species)
#'  * P_count: number of mosquitoes in the pupal stage (per species)
#'  * Sm_count: number of adult female mosquitoes who are Susceptible (per
#'  species)
#'  * Pm_count: number of adult female mosquitoes who are incubating (per
#'  species)
#'  * Im_count: number of adult female mosquitoes who are infectious (per
#'  species)
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
      correlations,
      list(create_lagged_infectivity(variables, parameters))
    ),
    variables = variables,
    events = unlist(events),
    timesteps = timesteps
  )
  renderer$to_dataframe()
}

#' @title Run a metapopulation model
#'
#' @param timesteps the number of timesteps to run the simulation for (in days)
#' @param parameters a named list of parameters to use
#' @param correlations correlation parameters
#' @param mixing matrix of mixing coefficients for infectivity
#' @return dataframe of results
#' @export
run_metapop_simulation <- function(
  timesteps,
  parameters,
  correlations = NULL,
  mixing
  ) {
  random_seed(ceiling(runif(1) * .Machine$integer.max))
  if (is.null(correlations)) {
    correlations <- lapply(parameters, get_correlation_parameters)
  }
  variables <- lapply(parameters, create_variables)
  events <- lapply(parameters, create_events)
  renderer <- lapply(parameters, function(.) individual::Render$new(timesteps))
  for (i in seq_along(parameters)) {
    initialise_events(events[[i]], variables[[i]], parameters[[i]])
    attach_event_listeners(
      events[[i]],
      variables[[i]],
      parameters[[i]],
      correlations[[i]],
      renderer[[i]]
    )
  }
  vector_models <- lapply(parameters, parameterise_mosquito_models)
  solvers <- lapply(
    seq_along(parameters),
    function(i) parameterise_solvers(vector_models[[i]], parameters[[i]])
  )
  lagged_infectivity <- lapply(
    seq_along(parameters),
    function(i) create_lagged_infectivity(variables[[i]], parameters[[i]])
  )
  processes <- lapply(
    seq_along(parameters),
    function(i) {
      create_processes(
        renderer[[i]],
        variables[[i]],
        events[[i]],
        parameters[[i]],
        vector_models[[i]],
        solvers[[i]],
        correlations[[i]],
        lagged_infectivity,
        mixing[i,],
        i
      )
    }
  )
  individual::simulation_loop(
    processes = unlist(processes),
    variables = unlist(variables),
    events = unlist(events),
    timesteps = timesteps
  )
  
  lapply(renderer, function(r) r$to_dataframe())
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

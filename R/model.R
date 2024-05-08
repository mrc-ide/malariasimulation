#' @title Run the simulation
#'
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
#'  * Tr_count: number of detectable infections being treated in humans
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
#'  * net_usage: the number people protected by a bed net
#'  * mosquito_deaths: number of adult female mosquitoes who die this timestep
#'  * n_drug_efficacy_failures: number of clinically treated individuals whose treatment failed due to drug efficacy
#'  * n_early_treatment_failure: number of clinically treated individuals who experienced early treatment failure
#'  * n_successfully_treated: number of clinically treated individuals who are treated successfully (includes individuals who experience slow parasite clearance)
#'  * n_slow_parasite_clearance: number of clinically treated individuals who experienced slow parasite clearance
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
  run_resumable_simulation(timesteps, parameters, correlations)$data
}

#' @title Run the simulation in a resumable way
#'
#' @description this function accepts an initial simulation state as an argument, and returns the
#' final state after running all of its timesteps. This allows one run to be resumed, possibly
#' having changed some of the parameters.
#' @param timesteps the timestep at which to stop the simulation
#' @param parameters a named list of parameters to use
#' @param correlations correlation parameters
#' @param initial_state the state from which the simulation is resumed
#' @param restore_random_state if TRUE, restore the random number generator's state from the checkpoint.
#' @return a list with two entries, one for the dataframe of results and one for the final
#' simulation state.
run_resumable_simulation <- function(
    timesteps,
    parameters = NULL,
    correlations = NULL,
    initial_state = NULL,
    restore_random_state = FALSE
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
  vector_models <- parameterise_mosquito_models(parameters, timesteps)
  solvers <- parameterise_solvers(vector_models, parameters)

  lagged_eir <- create_lagged_eir(variables, solvers, parameters)
  lagged_infectivity <- create_lagged_infectivity(variables, parameters)

  stateful_objects <- list(
    RandomState$new(restore_random_state),
    correlations,
    vector_models,
    solvers,
    lagged_eir,
    lagged_infectivity)

  if (!is.null(initial_state)) {
    individual::restore_object_state(
      initial_state$timesteps,
      stateful_objects,
      initial_state$malariasimulation)
  }

  individual_state <- individual::simulation_loop(
    processes = create_processes(
      renderer,
      variables,
      events,
      parameters,
      vector_models,
      solvers,
      correlations,
      list(lagged_eir),
      list(lagged_infectivity),
      timesteps
    ),
    variables = variables,
    events = events,
    timesteps = timesteps,
    state = initial_state$individual,
    restore_random_state = restore_random_state
  )

  final_state <- list(
    timesteps = timesteps,
    individual = individual_state,
    malariasimulation = individual::save_object_state(stateful_objects)
  )

  data <- renderer$to_dataframe()
  if (!is.null(initial_state)) {
    # Drop the timesteps we didn't simulate from the data.
    # It would just be full of NA.
    data <- data[-(1:initial_state$timesteps),]
  }

  list(data=data, state=final_state)
}

#' @title Run a metapopulation model
#'
#' @param timesteps the number of timesteps to run the simulation for (in days)
#' @param parameters a list of model parameter lists for each population
#' @param correlations a list of correlation parameters for each population
#' (default: NULL)
#' @param mixing matrix of mixing coefficients for infectivity towards
#' mosquitoes. Rows = origin sites, columns = destinations. Each element must 
#' be between 0 and 1 and all rows must sum to 1.
#' @return a list of dataframe of results
#' @export
run_metapop_simulation <- function(
    timesteps,
    parameters,
    correlations = NULL,
    mixing
) {
  random_seed(ceiling(runif(1) * .Machine$integer.max))
  if (nrow(mixing) != ncol(mixing)) {
    stop('mixing matrix must be square')
  }
  if (nrow(mixing) != length(parameters)) {
    stop('mixing matrix rows must match length of parameters')
  }
  if (!all(vlapply(seq_along(parameters), function(x) approx_sum(mixing[x,], 1)))) {
    stop('all mixing matrix rows must sum to 1')
  }
  if (!all(vlapply(seq_along(parameters), function(x) approx_sum(mixing[,x], 1)))) {
    warning('mixing matrix is asymmetrical')
  }
  if (is.null(correlations)) {
    correlations <- lapply(parameters, get_correlation_parameters)
  }
  variables <- lapply(parameters, create_variables)
  events <- lapply(parameters, create_events)
  renderer <- lapply(parameters, function(.) individual::Render$new(timesteps))
  for (i in seq_along(parameters)) {
    # NOTE: forceAndCall is necessary here to make sure i refers to the current
    # iteration
    forceAndCall(
      3,
      initialise_events,
      events[[i]],
      variables[[i]],
      parameters[[i]]
    )
    forceAndCall(
      5,
      attach_event_listeners,
      events[[i]],
      variables[[i]],
      parameters[[i]],
      correlations[[i]],
      renderer[[i]]
    )
  }
  vector_models <- lapply(parameters, parameterise_mosquito_models, timesteps = timesteps)
  solvers <- lapply(
    seq_along(parameters),
    function(i) parameterise_solvers(vector_models[[i]], parameters[[i]])
  )
  lagged_eir <- lapply(
    seq_along(parameters),
    function(i) create_lagged_eir(variables[[i]], solvers[[i]], parameters[[i]])
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
        lagged_eir,
        lagged_infectivity,
        timesteps,
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

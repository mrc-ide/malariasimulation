ageing_process <- function(simulation_frame, timestep, parameters) {
  if (timestep %% (365*parameters$timestep_to_day) == 0) {
    return(
      VariableUpdate$new(
        human,
        age,
        simulation_frame$get_variable(human, age) + 1
      )
    )
  }
}

maternal_immunity_process <- function(simulation_frame, timestep, parameters) {
  i <- simulation_frame$get_variable(human, icm)
  VariableUpdate$new(human, icm, i - parameters$rm * i)
}

preerythoctic_immunity_process <- function(simulation_frame, timestep, parameters) {
  i <- simulation_frame$get_variable(human, ib)
  VariableUpdate$new(
    human,
    ib,
    immunity_decay(
      i,
      simulation_frame$get_variable(human, last_bitten),
      timestep,
      parameters$rb,
      parameters$ub
    )
  )
}

acquired_immunity_process <- function(simulation_frame, timestep, parameters) {
  i <- simulation_frame$get_variable(human, ica)
  VariableUpdate$new(
    human,
    ica,
    immunity_decay(
      i,
      simulation_frame$get_variable(human, last_infected),
      timestep,
      parameters$rc,
      parameters$uc
    )
  )
}

infection_process <- function(simulation_frame, timestep, parameters) {
  source_humans <- simulation_frame$get_state(human, S, U, A)
  source_mosquitos <- simulation_frame$get_state(mosquito, Im)

  lambda <- force_of_infection(
    simulation_frame$get_variable(human, age)[source_humans],
    simulation_frame$get_constant(human, xi)[source_humans],
    simulation_frame$get_constant(mosquito, mosquito_variety)[source_mosquitos],
    simulation_frame$get_variable(human, ib)[source_humans],
    parameters
  )

  infected_humans <- source_humans[runif(length(source_humans), 0, 1) > lambda]

  phi <- immunity(
    simulation_frame$get_variable(human, ica)[infected_humans],
    simulation_frame$get_variable(human, icm)[infected_humans],
    parameters
  )

  symptomatic <- runif(length(infected_humans), 0, 1) > phi

  list(
    StateUpdate$new(human, I, infected_humans[symptomatic]),
    StateUpdate$new(human, A, infected_humans[!symptomatic]),
    VariableUpdate$new(
      human,
      last_bitten,
      timestep,
      infected_humans
    ),
    VariableUpdate$new(
      human,
      last_infected,
      timestep,
      infected_humans[symptomatic]
    )
  )
}

mosquito_infection_process <- function(simulation_frame, timestep, parameters) {
  source_mosquitos <- simulation_frame$get_state(mosquito, Sm)
  human_age <- simulation_frame$get_variable(human, age)
  human_xi <- simulation_frame$get_constant(human, xi)

  # Create a dataframe frame with age, xi and infectivity
  human_frame <- do.call(
    'rbind',
    lapply(
      list(
        list(D, parameters$cd),
        list(Treated, parameters$ct),
        list(A, parameters$ca),
        list(U, parameters$cu)
      ),
      function(args) {
        subset <- simulation_frame$get_state(human, args[[1]])
        if (length(subset) == 0) {
          return(
            setNames(
              data.frame(
                matrix(ncol = 3, nrow = 0)
              ),
              c('age', 'xi', 'infectivity')
            )
          )
        }
        data.frame(
          age         = human_age[subset],
          xi          = human_xi[subset],
          infectivity = args[[2]]
        )
      }
    )
  )

  lambda <- mosquito_force_of_infection(
    simulation_frame$get_constant(mosquito, mosquito_variety)[source_mosquitos],
    human_frame,
    parameters
  )
  infected = source_mosquitos[
    runif(length(source_mosquitos), 0, 1) > lambda
  ]
  StateUpdate$new(mosquito, Im, infected)
}

bind_process_to_model <- function(process, individuals, states, variables) {
  env <- environment(process)
  list2env(individuals, env)
  list2env(states, env)
  list2env(variables, env)
  process
}

create_processes <- function(individuals, states, variables, parameters) {
  
  env <- environment()
  list2env(individuals, env)
  list2env(states, env)

  processes <- list(
    # ===============
    # Human Processes
    # ===============

    ageing_process,
    maternal_immunity_process,
    preerythoctic_immunity_process,
    acquired_immunity_process,

    # Untreated Progression
    create_fixed_probability_state_change_process(human, I, D, 1 - parameters$ft),
    # Treatment
    create_fixed_probability_state_change_process(human, I, Treated, parameters$ft),
    # Asymptomatic Progression
    create_fixed_probability_state_change_process(human, D, A, parameters$rd),
    # Subpatient Progression
    create_fixed_probability_state_change_process(human, A, U, parameters$ra),
    # Subpatient Recovery
    create_fixed_probability_state_change_process(human, U, S, parameters$ru),
    # Treatment Recovery
    create_fixed_probability_state_change_process(human, Treated, S, parameters$rt),

    # ==================
    # Mosquito Processes
    # ==================

    # Larval growth
    create_fixed_probability_state_change_process(mosquito, E, L, parameters$rel),
    # Pupal stage
    create_fixed_probability_state_change_process(mosquito, L, P, parameters$rl),
    # Susceptable Female Development
    create_fixed_probability_state_change_process(mosquito, P, Sm, parameters$rpl),
    ## Mosquito Infection
    ## Mosquitos move from Sm -> Im
    ## NOTE: In the future this will be combined with the infection process
    ## below so that we can model mosquitos individually
    mosquito_infection_process,

    ## =============
    ## Mosquito bite
    ## =============
    ## Mosquito bites move humans from S, U, A -> I; S, U, A -> A
    ## and has a side effect of boosting immunity
    infection_process
  )

  lapply(
    processes,
    function(p) { bind_process_to_model(p, individuals, states, variables) }
  )
}

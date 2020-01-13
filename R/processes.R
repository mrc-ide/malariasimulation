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

create_processes <- function(individuals, states, variables, parameters) {
  
  env <- environment()
  list2env(individuals, env)
  list2env(states, env)
  list2env(variables, env)

  processes <- list(
    ageing_process,

    # ======
    # States
    # ======
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

    # ========
    # Immunity
    # ========
    # Maternal immunity
    create_exponential_decay_process(human, icm, parameters$rm),
    create_exponential_decay_process(human, ivm, parameters$rvm),
    # Blood immunity
    create_exponential_decay_process(human, ib, parameters$rb),
    # Acquired immunity
    create_exponential_decay_process(human, ica, parameters$rc),
    create_exponential_decay_process(human, iva, parameters$rva),

    # ===============
    # Mosquito States
    # ===============

    # Larval growth
    create_fixed_probability_state_change_process(mosquito, E, L, parameters$rel),
    # Pupal stage
    create_fixed_probability_state_change_process(mosquito, L, P, parameters$rl),
    # Susceptable Female Development
    create_fixed_probability_state_change_process(mosquito, P, Sm, parameters$rpl),

    ## =========
    ## Infection
    ## =========
    ## Mosquitos move from Sm -> Im
    ## NOTE: In the future this will be combined with the infection process
    ## below so that we can model mosquitos individually
    mosquito_infection_process,

    ## schedule infections for humans and set last_bitten and last_infected
    infection_process,

    ## update states after a latent period
    scheduled_infections
  )

  lapply(
    processes,
    function(p) { bind_process_to_model(p, individuals, states, variables) }
  )
}

# =================
# Utility functions
# =================

bind_process_to_model <- function(process, individuals, states, variables) {
  env <- environment(process)
  list2env(individuals, env)
  list2env(states, env)
  list2env(variables, env)
  process
}

create_fixed_probability_state_change_process <- function(i, from, to, rate) {
  function (simulation_frame, timestep, parameters) {
    source_individuals <- simulation_frame$get_state(i, from)
    target_individuals <- source_individuals[
      runif(length(source_individuals), 0, 1) > rate
    ]
    StateUpdate$new(i, to, target_individuals)
  }
}

create_exponential_decay_process <- function(individual, variable, rate) {
  function(simulation_frame, timestep, parameters) {
    i <- simulation_frame$get_variable(individual, variable)
    VariableUpdate$new(individual, variable, i - rate * i)
  }
}

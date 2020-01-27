
#' create_processes, defines the functions which describe how each individual's
#' states and variables change over time.
#'
#' It lists processes from `infection.R`, `mosquito_emergence.R` and
#' `mortality.R`; and then exposes them to the model through
#' `bind_process_to_model`
#' @param individuals, a list of individuals in the model
#' @param states, a list of states in the model
#' @param variables, a list of variables and constants in the model
#' @param parameters, a list of model parameters
create_processes <- function(individuals, states, variables, parameters) {
  
  env <- environment()
  list2env(individuals, env)
  list2env(states, env)
  list2env(variables, env)

  processes <- list(
    ageing_process,
    mortality_process,

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
    treatment_recovery_process,

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
    create_exponential_decay_process(human, id, parameters$rd),

    # ===============
    # Mosquito States
    # ===============
    # Eggs laid
    egg_laying_process,
    # Larval growth
    create_fixed_probability_state_change_process(mosquito, E, L, parameters$rel),
    # Pupal stage
    create_fixed_probability_state_change_process(mosquito, L, P, parameters$rl),
    # Susceptable Female Development
    create_fixed_probability_state_change_process(mosquito, P, Sm, .5 * parameters$rpl),
    # Death of larvae
    larval_death_process,
    # Death of pupals
    create_fixed_probability_state_change_process(mosquito, P, Unborn, .5 * parameters$mup),
    # Natural death of females
    create_fixed_probability_state_change_process(mosquito, Sm, Unborn, parameters$mum),
    create_fixed_probability_state_change_process(mosquito, Im, Unborn, parameters$mum),

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

#' bind_process_to_model adds individuals, states and variables to a process
#' functions's environment so that it can specify model updates at each timestep
#'
#' @param individuals, a list of individuals in the model
#' @param states, a list of states in the model
#' @param variables, a list of variables and constants in the model
bind_process_to_model <- function(process, individuals, states, variables) {
  env <- environment(process)
  list2env(individuals, env)
  list2env(states, env)
  list2env(variables, env)
  process
}

#' create_fixed_probability_state_change_process generates a process function
#' that moves individuals from one state to another at a constant rate
#'
#' @param i, an individual
#' @param from, the source state
#' @param to, the target state
#' @param rate, the rate at which state transitions occur
create_fixed_probability_state_change_process <- function(i, from, to, rate) {
  function (simulation_frame, timestep, parameters) {
    source_individuals <- simulation_frame$get_state(i, from)
    target_individuals <- source_individuals[
      bernoulli(length(source_individuals), rate)
    ]
    individual::StateUpdate$new(i, to, target_individuals)
  }
}

#' create_exponential_decay_process generates a process function
#' that reduces the value of a variable at an exponential rate
#'
#' @param i, an individual
#' @param variable, the variable to update
#' @param rate, the exponential rate
create_exponential_decay_process <- function(individual, variable, rate) {
  function(simulation_frame, timestep, parameters) {
    i <- simulation_frame$get_variable(individual, variable)
    individual::VariableUpdate$new(individual, variable, i - rate * i)
  }
}

#' This is the process for aging, it will update every human's age every 365
#' timesteps.
#'
#' @param simulation_frame, the current state of the simulation
#' @param timestep, the current timestep
#' @param parameters, the model parameters
ageing_process <- function(simulation_frame, timestep, parameters) {
  if (timestep %% (365*parameters$timestep_to_day) == 0) {
    return(
      individual::VariableUpdate$new(
        human,
        age,
        simulation_frame$get_variable(human, age) + 1
      )
    )
  }
}

#' @title Define model processes
#' @description
#' create_processes, defines the functions which describe how each individual's
#' states and variables change over time.
#'
#' It lists processes from `infection.R`, `mosquito_emergence.R` and
#' `mortality.R`; and then exposes them to the model through
#' `bind_process_to_model`
#' @param individuals, a list of individuals in the model
#' @param states, a list of states in the model
#' @param variables, a list of variables in the model
#' @param parameters, a list of model parameters
create_processes <- function(individuals, states, variables, parameters) {
  list(
    create_ageing_process(individuals$human, variables$age),
    create_mortality_process(
      individuals$human,
      states$D,
      variables
    ),

    # ======
    # States
    # ======

    # Untreated Progression
    create_fixed_probability_state_change_process(
      individuals$human,
      states$I,
      states$D,
      1 - parameters$ft
    ),
    # Asymptomatic Progression
    create_fixed_probability_state_change_process(
      individuals$human,
      states$D,
      states$A,
      parameters$rd
    ),
    # Subpatient Progression
    create_fixed_probability_state_change_process(
      individuals$human,
      states$A,
      states$U,
      parameters$ra
    ),
    # Subpatient Recovery
    create_fixed_probability_state_change_process(
      individuals$human,
      states$U,
      states$S,
      parameters$ru
    ),

    # ========
    # Immunity
    # ========

    # Maternal immunity
    create_exponential_decay_process(individuals$human, variables$icm, parameters$rm),
    create_exponential_decay_process(individuals$human, variables$ivm, parameters$rvm),
    # Blood immunity
    create_exponential_decay_process(individuals$human, variables$ib, parameters$rb),
    # Acquired immunity
    create_exponential_decay_process(individuals$human, variables$ica, parameters$rc),
    create_exponential_decay_process(individuals$human, variables$iva, parameters$rva),
    create_exponential_decay_process(individuals$human, variables$id, parameters$rd),

    # ========
    # Immunity
    # ========

    # Maternal immunity
    create_exponential_decay_process(individuals$human, variables$icm, parameters$rm),
    create_exponential_decay_process(individuals$human, variables$ivm, parameters$rvm),
    # Blood immunity
    create_exponential_decay_process(individuals$human, variables$ib, parameters$rb),
    # Acquired immunity
    create_exponential_decay_process(individuals$human, variables$ica, parameters$rc),
    create_exponential_decay_process(individuals$human, variables$iva, parameters$rva),
    create_exponential_decay_process(individuals$human, variables$id, parameters$rid),

    # ===============
    # Mosquito States
    # ===============
    # Eggs laid
    create_egg_laying_process(
      individuals$mosquito,
      states$Sm,
      states$Im,
      states$Unborn,
      states$E
    ),
    # Larval growth
    create_fixed_probability_state_change_process(
      individuals$mosquito,
      states$E,
      states$L,
      parameters$rel
    ),
    # Pupal stage
    create_fixed_probability_state_change_process(
      individuals$mosquito,
      states$L,
      states$P,
      parameters$rl
    ),
    # Susceptable Female Development
    create_fixed_probability_state_change_process(
      individuals$mosquito,
      states$P,
      states$Sm,
      .5 * parameters$rpl
    ),
    # Death of larvae
    create_larval_death_process(
      individuals$mosquito,
      states$E,
      states$L,
      states$Unborn
    ),
    # Death of pupals
    create_fixed_probability_state_change_process(
      individuals$mosquito,
      states$P,
      states$Unborn,
      parameters$mup
    ),
    # Natural death of females
    create_fixed_probability_state_change_process(
      individuals$mosquito,
      states$Sm,
      states$Unborn,
      parameters$mum
    ),
    create_fixed_probability_state_change_process(
      individuals$mosquito,
      states$Im,
      states$Unborn,
      parameters$mum
    ),

    # =========
    # Infection
    # =========
    create_mosquito_infection_process(
      individuals$mosquito,
      individuals$human,
      states,
      variables
    ),

    # schedule infections for humans and set last_bitten and last_infected
    create_infection_process(
      individuals$human,
      individuals$mosquito,
      states,
      variables
    ),

    # update states after a latent period
    create_infection_scheduler(
      individuals$human,
      states$A,
      states$I,
      variables$infection_schedule,
      variables$asymptomatic_infection_schedule
    )
  )
}

# =================
# Utility functions
# =================

#' @title Basic state transition
#' @description
#' create_fixed_probability_state_change_process generates a process function
#' that moves individuals from one state to another at a constant rate
#'
#' @param i, an individual
#' @param from, the source state
#' @param to, the target state
#' @param rate, the rate at which state transitions occur
create_fixed_probability_state_change_process <- function(i, from, to, rate) {
  stopifnot(is.numeric(rate))
  function (simulation_frame, timestep, parameters) {
    source_individuals <- simulation_frame$get_state(i, from)
    target_individuals <- source_individuals[
      bernoulli(length(source_individuals), rate)
    ]
    individual::StateUpdate$new(i, to, target_individuals)
  }
}

#' @title Exponentially decaying variables
#' @description
#' create_exponential_decay_process generates a process function
#' that reduces the value of a variable at an exponential rate
#'
#' @param individual, an individual
#' @param variable, the variable to update
#' @param rate, the exponential rate
create_exponential_decay_process <- function(individual, variable, rate) {
  function(simulation_frame, timestep, parameters) {
    i <- simulation_frame$get_variable(individual, variable)
    individual::VariableUpdate$new(individual, variable, i - rate * i)
  }
}

#' @title Human aging process
#' @description
#' This is the process for aging, it will update every human's age every 365
#' timesteps.
#'
#' @param human, the human individual
#' @param age, the age variable
create_ageing_process <- function(human, age) {
  function(simulation_frame, timestep, parameters) {
    if (timestep %% (365 / parameters$days_per_timestep) == 0) {
      return(
        individual::VariableUpdate$new(
          human,
          age,
          simulation_frame$get_variable(human, age)[0] + 1
        )
      )
    }
  }
}

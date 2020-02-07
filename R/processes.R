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
#' @param variables, a list of variables and constants in the model
#' @param parameters, a list of model parameters
create_processes <- function(individuals, states, variables, parameters) {

  list(
    create_ageing_process(individuals$human, variables$age),

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
    # Treatment
    create_fixed_probability_state_change_process(
      individuals$human,
      states$I,
      states$Treated,
      parameters$ft
    ),
    # Prophylaxis
    create_fixed_probability_state_change_process(
      individuals$human,
      states$I,
      states$Treated,
      parameters$rp
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
    # Death of pupals
    create_fixed_probability_state_change_process(
      individuals$mosquito,
      states$P,
      states$Unborn,
      .5 * parameters$mup
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

#' @title Human aging process
#' @description
#' This is the process for aging, it will update every human's age every 365
#' timesteps.
#'
#' @param human, the human individual
#' @param age, the age variable
create_ageing_process <- function(human, age) {
  function(simulation_frame, timestep, parameters) {
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
}

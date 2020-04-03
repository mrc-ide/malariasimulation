#' @title Define model processes
#' @description
#' create_processes, defines the functions which describe how each individual's
#' states and variables change over time.
#'
#' It lists processes from `infection.R`, `mosquito_emergence.R` and
#' `mortality.R`; and then exposes them to the model
#' @param individuals, a list of individuals in the model
#' @param states, a list of states in the model
#' @param variables, a list of variables in the model
#' @param events, a list of events in the model
#' @param parameters, a list of model parameters
create_processes <- function(individuals, states, variables, events, parameters) {
  list(
    create_ageing_process(individuals$human, variables$age),
    create_mortality_process(
      individuals$human,
      states$D,
      variables,
      events
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
    create_exponential_decay_process(individuals$human, variables$id, parameters$rid),

    # schedule infections for humans and set last_bitten and last_infected
    create_infection_process(
      individuals,
      states,
      variables,
      events
    ),

    # ===============
    # Mosquito States
    # ===============
    create_mosquito_infection_process(
      individuals$mosquito,
      individuals$human,
      states,
      variables
    ),

    # Eggs laid
    create_egg_laying_process(
      individuals$mosquito,
      states,
      events,
      parameters
    ),

    # Death of larvae
    create_larval_death_process(
      individuals$mosquito,
      states$E,
      states$L,
      states$Unborn,
      events
    ),
    # Death of pupals
    create_pupal_death_process(
      individuals$mosquito,
      states$P,
      states$Unborn,
      parameters$mup,
      events
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

#' @title Define event based processes
#' @description
#' defines processes for events that can be scheduled in the future
#'
#' @param individuals, a list of individuals in the model
#' @param states, a list of states in the model
#' @param events, a list of events in the model
create_event_based_processes <- function(individuals, states, events, parameters) {

  # Disease progression events
  events$infection$add_listener(function(api, target) {
    api$schedule(events$asymptomatic_progression, target, parameters$dd)
    individual::StateUpdate$new(individuals$human, states$D, target)
  })
  events$asymptomatic_infection$add_listener(function(api, target) {
    api$schedule(events$subpatent_progression, target, parameters$da)
    individual::StateUpdate$new(individuals$human, states$A, target)
  })
  events$subpatent_progression$add_listener(function(api, target) {
    api$schedule(events$subpatent_recovery, target, parameters$du)
    individual::StateUpdate$new(individuals$human, states$U, target)
  })
  events$subpatent_recovery$add_listener(function(api, target) {
    individual::StateUpdate$new(individuals$human, states$S, target)
  })

  # Mosquito development processes
  events$larval_growth$add_listener(function(api, target) {
    api$schedule(events$pupal_development, target, parameters$dl)
    individual::StateUpdate$new(individuals$mosquito, states$L, target)
  })
  events$pupal_development$add_listener(function(api, target) {
    api$schedule(events$susceptable_development, target, parameters$dpl)
    individual::StateUpdate$new(individuals$mosquito, states$P, target)
  })
  events$susceptable_development$add_listener(function(api, target) {
    female <- bernoulli(length(target), .5)
    list(
      individual::StateUpdate$new(individuals$mosquito, states$Unborn, target[!female]),
      individual::StateUpdate$new(individuals$mosquito, states$Sm, target[female])
    )
  })
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
  function (api) {
    source_individuals <- api$get_state(i, from)
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
  function(api) {
    i <- api$get_variable(individual, variable)
    individual::VariableUpdate$new(individual, variable, i - rate * i)
  }
}

#TODO: Schedule the aging process so that it's more granular

#' @title Human aging process
#' @description
#' This is the process for aging, it will update every human's age every 365
#' timesteps.
#'
#' @param human, the human individual
#' @param age, the age variable
create_ageing_process <- function(human, age) {
  function(api) {
    if (api$get_timestep() %% (365 / api$get_parameters()$days_per_timestep) == 0) {
      return(
        individual::VariableUpdate$new(
          human,
          age,
          api$get_variable(human, age) + 1
        )
      )
    }
  }
}

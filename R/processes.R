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
    create_egg_laying_process_cpp(
      individuals$mosquito$name,
      states$Sm$name,
      states$Im$name,
      states$Unborn$name,
      states$E$name,
      events$larval_growth$name
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
    individual::fixed_probability_state_change_process(
      individuals$mosquito$name,
      states$Sm$name,
      states$Unborn$name,
      parameters$mum
    ),
    individual::fixed_probability_state_change_process(
      individuals$mosquito$name,
      states$Im$name,
      states$Unborn$name,
      parameters$mum
    ),

    # Rendering processes
    individual::state_count_renderer_process(
      individuals$human$name,
      c(states$S$name, states$A$name, states$D$name, states$U$name)
    ),
    individual::state_count_renderer_process(
      individuals$mosquito$name,
      c(states$E$name, states$L$name, states$P$name, states$Sm$name, states$Im$name)
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
create_event_based_processes <- function(individuals, states, variables, events, parameters) {
  # Aging
  events$birthday$add_listener(function(api, target) {
    api$schedule(events$birthday, target, 365)
    api$queue_variable_update(
      individuals$human,
      variables$age,
      api$get_variable(individuals$human, variables$age)[target] + 1,
      target
    )
  })

  # Disease progression events
  events$infection$add_listener(
    individual::update_state_listener(individuals$human$name, states$D$name)
  )
  events$infection$add_listener(
    individual::reschedule_listener(events$asymptomatic_progression$name, parameters$dd)
  )
  events$asymptomatic_infection$add_listener(
    individual::update_state_listener(individuals$human$name, states$A$name)
  )
  events$asymptomatic_infection$add_listener(
    individual::reschedule_listener(events$subpatent_progression$name, parameters$da)
  )
  events$subpatent_progression$add_listener(
    individual::update_state_listener(individuals$human$name, states$U$name)
  )
  events$subpatent_progression$add_listener(
    individual::reschedule_listener(events$subpatent_recovery$name, parameters$du)
  )
  events$subpatent_recovery$add_listener(
    individual::update_state_listener(individuals$human$name, states$S$name)
  )

  # Mosquito development processes
  events$larval_growth$add_listener(
    individual::update_state_listener(individuals$mosquito$name, states$L$name)
  )
  events$larval_growth$add_listener(
    individual::reschedule_listener(events$pupal_development$name, parameters$dl)
  )
  events$pupal_development$add_listener(
    individual::update_state_listener(individuals$mosquito$name, states$P$name)
  )
  events$pupal_development$add_listener(
    individual::reschedule_listener(events$susceptable_development$name, parameters$dpl)
  )
  events$susceptable_development$add_listener(function(api, target) {
    female <- bernoulli(length(target), .5)
    api$queue_state_update(individuals$mosquito, states$Unborn, target[!female])
    api$queue_state_update(individuals$mosquito, states$Sm, target[female])
  })
}

# =================
# Utility functions
# =================

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
    api$queue_variable_update(individual, variable, i - rate * i)
  }
}

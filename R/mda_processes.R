
#' @title Create listeners for MDA events
#' @param human the human individual object
#' @param states the states available in the model
#' @param variables the variables available in the model
#' @param mda_index the index of the MDA to create listeners for
#' @param mda_events a list of events for the MDA at `mda_index`
#' @description will create an "enrollment listener" for setting up a target
#' population for an MDA and an  "administer listener" for modelling the effects
#' of the drug on the population state
create_mda_listeners <- function(human, states, variables, mda_index, mda_events) {
  mda_administer_listener <- function(api, target) {
    parameters <- api$get_parameters()
    successful_treatments <- bernoulli(
      length(target),
      parameters$drug_efficacies[[
        parameters$mda_drug[[mda_index]]
      ]]
    )
    to_move <- target[successful_treatments]

    # Move Diseased
    diseased <- api$get_state(human, states$D, states$A)
    api$queue_state_update(
      human,
      states$T,
      intersect(to_move, diseased)
    )

    # Move everyone else
    api$queue_state_update(
      human,
      states$Ph,
      setdiff(to_move, diseased)
    )

    # Schedule next dose
    frequency <- parameters$mda_frequency[[mda_index]]
    if (api$get_timestep() + frequency <= parameters$mda_end[[mda_index]]) {
      api$schedule(
        mda_events$mda_administer,
        target,
        frequency
      )
    }
  }

  mda_enrollment_listener <- function(api, target) {
    parameters <- api$get_parameters()
    age <- get_age(
      api$get_variable(human, variables$birth),
      api$get_timestep()
    )
    target <- which(
      (age > parameters$mda_min_age[[mda_index]]) &
      (age < parameters$mda_max_age[[mda_index]])
    )
    covered <- bernoulli(length(target), parameters$mda_coverage[[mda_index]])
    mda_administer_listener(api, target[covered], mda_index)
  }

  list(
    mda_enrollment_listener = mda_enrollment_listener,
    mda_administer_listener = mda_administer_listener
  )
}

attach_mda_listeners <- function(human, states, variables, mda_events) {
  for (mda_index in seq_along(mda_events)) {
    listeners <- create_mda_listeners(
      human,
      states,
      variables,
      mda_index,
      mda_events[[mda_index]]
    )
    mda_events[[mda_index]]$mda_enrollment$add_listener(
      listeners$mda_enrollment_listener
    )
    mda_events[[mda_index]]$mda_administer$add_listener(
      listeners$mda_administer_listener
    )
  }
}

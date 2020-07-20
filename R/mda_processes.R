
#' @title Create listeners for MDA events
#' @param human the human individual object
#' @param states the states available in the model
#' @param variables the variables available in the model
#' @param administer_event the event schedule for drug administration
#' @param drug the drug to administer
#' @param end when to stop distributing (timestep)
#' @param frequency how often doses are distributed (in timesteps)
#' @param min_age minimum age at enrollment
#' @param max_age maximum age at enrollment
#' @param coverage the proportion of the target population that is covered
#' @description will create an "enrollment listener" for setting up a target
#' population for an MDA and an  "administer listener" for modelling the effects
#' of the drug on the population state
create_mda_listeners <- function(
  human,
  states,
  variables,
  administer_event,
  drug,
  end,
  frequency,
  min_age,
  max_age,
  coverage
  ) {

  administer_listener <- function(api, target) {
    parameters <- api$get_parameters()
    timestep <- api$get_timestep()
    successful_treatments <- bernoulli(
      length(target),
      parameters$drug_efficacies[[drug]]
    )
    to_move <- target[successful_treatments]

    if (length(to_move > 0)) {
      # Move Diseased
      diseased <- intersect(to_move, api$get_state(human, states$D, states$A))
      if (length(diseased) > 0) {
        api$queue_state_update(
          human,
          states$Tr,
          diseased
        )
      }

      # Move everyone else
      other <- setdiff(to_move, diseased)
      if (length(other) > 0) {
        api$queue_state_update(
          human,
          states$S,
          other
        )
      }

      # Update infectivity
      api$queue_variable_update(
        human,
        variables$infectivity,
        api$get_variable(
          human,
          variables$infectivity,
          to_move
        ) * parameters$drug_rel_c[[drug]],
        to_move
      )

      # Update drug
      api$queue_variable_update(human, variables$drug, drug, to_move)
      api$queue_variable_update(human, variables$drug_time, timestep, to_move)
    }

    # Schedule next dose
    if (timestep + frequency <= end) {
      api$schedule(
        administer_event,
        target,
        frequency
      )
    }
  }

  enrollment_listener <- function(api, target) {
    parameters <- api$get_parameters()
    age <- get_age(
      api$get_variable(human, variables$birth),
      api$get_timestep()
    )
    target <- which((age > min_age) & (age < max_age))
    covered <- bernoulli(length(target), coverage)
    administer_listener(api, target[covered])
  }

  list(
    enrollment_listener = enrollment_listener,
    administer_listener = administer_listener
  )
}

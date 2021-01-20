#' @title Create listeners for MDA events
#' @param variables the variables available in the model
#' @param administer_event the event schedule for drug administration
#' @param drug the drug to administer
#' @param end when to stop distributing (timestep)
#' @param frequency how often doses are distributed (in timesteps)
#' @param min_age minimum age for the target population
#' @param max_age maximum age for the target population
#' @param coverage the proportion of the target population that is covered
#' @param correlations correlation parameters
#' @param int_name the name of this intervention (either 'smc' or 'mda')
#' @description will create a listener for administering each round of drugs
create_mda_listeners <- function(
  variables,
  administer_event,
  drug,
  end,
  frequency,
  min_age,
  max_age,
  coverage,
  correlations,
  int_name,
  parameters
  ) {
  function(timestep, target) {
    age <- get_age(variables$birth$get_values(), timestep)

    in_age <- which((age > min_age) & (age < max_age))
    target <- in_age[sample_intervention(in_age, int_name, coverage, correlations)]

    successful_treatments <- bernoulli(
      length(target),
      parameters$drug_efficacy[[drug]]
    )
    to_move <- target[successful_treatments]

    renderer$render('n_mda_treated', length(to_move))

    if (length(to_move > 0)) {
      # Move Diseased
      diseased <- intersect(to_move, variables$state$get_index_of('D', 'A')$to_vector())
      if (length(diseased) > 0) {
        variables$state$queue_update('Tr', diseased)
      }

      # Move everyone else
      other <- setdiff(to_move, diseased)
      if (length(other) > 0) {
        variables$state$queue_update('S', other)
      }

      # Update infectivity
      variables$infectivity$queue_update(
        variables$infectivity$get_values(
          human,
          to_move
        ) * parameters$drug_rel_c[[drug]],
        to_move
      )

      # Update drug
      variables$drug$queue_update(drug, to_move)
      variables$drug_time$queue_update(timestep, to_move)
    }

    # Schedule next round
    if (timestep + frequency <= end) {
      administer_event$schedule(frequency)
    }
  }
}

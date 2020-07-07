create_events <- function() {
  list(
    # Human infection events
    infection = individual::Event$new('infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection'),

    # Mosquito events
    mosquito_infection = individual::Event$new('mosquito_infection')
  )
}

create_mda_events <- function(parameters) {
  events <- list()

  for (mda_index in seq_along(parameters$mda_drug)) {
    events[[length(events) + 1]] <- list(
      mda_enrollment = individual::Event$new(paste0('mda_enrolment_', mda_index)),
      mda_administer = individual::Event$new(paste0('mda_administer_', mda_index))
    )
  }

  events
}

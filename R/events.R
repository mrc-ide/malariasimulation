create_events <- function() {
  list(
    # Human infection events
    infection = individual::Event$new('infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection'),

    # Human MDA events
    mda_enrollment = individual::Event$new('mda_enrolment'),
    mda_administer = individual::Event$new('mda_administer'),

    # Mosquito events
    mosquito_infection = individual::Event$new('mosquito_infection')
  )
}

create_events <- function() {
  list(
    # Human infection events
    infection = individual::Event$new('infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection'),

    # Vaccination events
    rtss_vaccination = individual::Event$new('rtss_vaccination'),

    # MDA events
    mda_enrollment = individual::Event$new('mda_enrolment'),
    mda_administer = individual::Event$new('mda_administer'),
    smc_enrollment = individual::Event$new('smc_enrolment'),
    smc_administer = individual::Event$new('smc_administer'),

    # Mosquito events
    mosquito_infection = individual::Event$new('mosquito_infection')
  )
}

create_events <- function() {
  list(
    # Human infection events
    clinical_infection = individual::Event$new('clinical_infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection'),
    infection = individual::Event$new('infection'),

    # Vaccination events
    rtss_vaccination = individual::Event$new('rtss_vaccination'),
    rtss_booster = individual::Event$new('rtss_booster'),

    # MDA events
    mda_enrollment = individual::Event$new('mda_enrolment'),
    mda_administer = individual::Event$new('mda_administer'),
    smc_enrollment = individual::Event$new('smc_enrolment'),
    smc_administer = individual::Event$new('smc_administer'),

    # Mosquito events
    mosquito_infection = individual::Event$new('mosquito_infection')
  )
}

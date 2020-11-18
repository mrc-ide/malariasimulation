create_events <- function() {
  list(
    # Human infection events
    clinical_infection = individual::Event$new('clinical_infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection'),
    infection = individual::Event$new('infection'), # either clinical or asym infection
    subpatent_infection = individual::Event$new('subpatent_infection'),
    recovery = individual::Event$new('recovery'),

    # Vaccination events
    rtss_vaccination = individual::Event$new('rtss_vaccination'),
    rtss_booster = individual::Event$new('rtss_booster'),
    tbv_vaccination = individual::Event$new('tbv_vaccination'),

    # MDA events
    mda_administer = individual::Event$new('mda_administer'),
    smc_administer = individual::Event$new('smc_administer'),

    # Bednet events
    throw_away_net = individual::Event$new('throw_away_net'),

    # Mosquito events
    mosquito_infection = individual::Event$new('mosquito_infection')
  )
}

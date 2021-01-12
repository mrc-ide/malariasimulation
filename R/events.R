create_events <- function(parameters) {
  list(
    # Human infection events
    clinical_infection = individual::TargetedEvent$new(parameters$human_population),
    asymptomatic_infection = individual::TargetedEvent$new(parameters$human_population),
    infection = individual::TargetedEvent$new(parameters$human_population), # either clinical or asym infection
    subpatent_infection = individual::TargetedEvent$new(parameters$human_population),
    recovery = individual::TargetedEvent$new(parameters$human_population),

    # Vaccination events
    rtss_vaccination = individual::Event$new(),
    rtss_booster = individual::TargetedEvent$new(parameters$human_population),
    tbv_vaccination = individual::Event$new(),

    # MDA events
    mda_administer = individual::Event$new(),
    smc_administer = individual::Event$new(),

    # Bednet events
    throw_away_net = individual::TargetedEvent$new(parameters$human_population),

    # Mosquito events
    mosquito_infection = individual::TargetedEvent$new(parameters$mosquito_limit)
  )
}

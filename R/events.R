create_events <- function() {
  list(
    # Human infection events
    infection = individual::Event$new('infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection'),

    # Vaccination events
    rtss_vaccination = individual::Event$new('rtss_vaccination'),

    mosquito_infection = individual::Event$new('mosquito_infection')
  )
}

create_events <- function() {
  list(
    # Human infection events
    infection = individual::Event$new('infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection'),

    mosquito_infection = individual::Event$new('mosquito_infection')
  )
}

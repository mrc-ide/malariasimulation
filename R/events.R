create_events <- function() {
  list(
    infection = individual::Event$new('infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection')
  )
}

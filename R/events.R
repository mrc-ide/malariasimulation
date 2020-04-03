create_events <- function() {
  list(
    # Human aging
    birthday = individual::Event$new('birthday'),

    # Human infection events
    infection = individual::Event$new('infection'),
    asymptomatic_infection = individual::Event$new('asymptomatic_infection'),
    asymptomatic_progression = individual::Event$new('asymptomatic_progression'),
    subpatent_progression = individual::Event$new('subpatent_progression'),
    subpatent_recovery = individual::Event$new('subpatent_recovery'),

    # Mosquito development events
    larval_growth = individual::Event$new('larval_growth'),
    pupal_development = individual::Event$new('pupal_development'),
    susceptable_development = individual::Event$new('susceptable_development')
  )
}

test_that('biting_process integrates mosquito effects and human infection', {
  population <- 4
  parameters <- get_parameters(list(human_population = population, severe_enabled = TRUE))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  biting_process <- create_biting_process(
    individuals,
    states,
    variables,
    events
  )

  age <- c(20, 24, 5, 39) * 365

  api <- mock_api(
    list(
      human = list(
        birth = -age + 5
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  total_eir <- 5
  eir <- rep(
    total_eir / population / length(parameters$blood_meal_rates),
    population
  )

  bites_mock <- mockery::mock(eir)
  infection_mock <- mockery::mock()

  mockery::stub(biting_process, 'simulate_bites', bites_mock)
  mockery::stub(biting_process, 'simulate_infection', infection_mock)
  biting_process(api)

  mockery::expect_args(
    bites_mock,
    1,
    api,
    individuals,
    states,
    variables,
    age,
    parameters
  )

  mockery::expect_args(
    infection_mock,
    1,
    api,
    individuals,
    states,
    variables,
    events,
    eir,
    age,
    parameters
  )
})

test_that('simulate_bites integrates eir calculation and mosquito side effects', {
  population <- 4
  parameters <- get_parameters(list(human_population = population, severe_enabled = TRUE))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  infectivity <- c(.6, 0, .2, .3)
  age <- c(20, 24, 5, 39) * 365

  api <- mock_api(
    list(
      human = list(
        zeta = c(.2, .3, .5, .9),
        infectivity = infectivity
      ),
      mosquito = list(
        Im = 1:10,
        Sm = 11:25,
        variety = c(rep(1, 25), rep(2, 25), rep(3, 50))
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  total_eir <- 5
  lambda_mock <- mockery::mock(.5, .5, .5)
  mosquito_effects_mock <- mockery::mock()

  mockery::stub(simulate_bites, 'effective_biting_rate', lambda_mock)
  mockery::stub(simulate_bites, 'calculate_mosquito_effects', mosquito_effects_mock)
  total_eir <- simulate_bites(api, individuals, states, variables, age, parameters)

  expect_equal(total_eir, 5)

  f <- parameters$blood_meal_rates[[1]]

  mockery::expect_args(
    lambda_mock,
    1,
    api,
    individuals$human,
    variables$zeta,
    age,
    1,
    list(
      prob_bitten_survives = rep(1, population),
      prob_bitten = rep(1, population),
      prob_repelled = rep(0, population)
    ),
    f,
    parameters
  )

  mockery::expect_args(
    mosquito_effects_mock,
    1,
    api,
    infectivity,
    .5,
    individuals,
    states,
    1,
    11:25,
    c(11:25, 1:10),
    1,
    0,
    f,
    parameters
  )
})

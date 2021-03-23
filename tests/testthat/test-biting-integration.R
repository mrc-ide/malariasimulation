test_that('biting_process integrates mosquito effects and human infection', {
  population <- 4
  timestep <- 5
  parameters <- get_parameters(
    list(human_population = population, severe_enabled = TRUE)
  )

  renderer <- individual::Render$new(5)
  events <- mockery::mock()
  age <- c(20, 24, 5, 39) * 365
  variables <- list(birth = individual::DoubleVariable$new((-age + timestep)))

  biting_process <- create_biting_process(
    renderer,
    list(), #solvers
    list(), #models
    variables,
    events,
    parameters
  )

  bitten <- individual::Bitset$new(parameters$human_population)
  bites_mock <- mockery::mock(bitten)
  infection_mock <- mockery::mock()

  mockery::stub(biting_process, 'simulate_bites', bites_mock)
  mockery::stub(biting_process, 'simulate_infection', infection_mock)
  biting_process(timestep)

  mockery::expect_args(
    bites_mock,
    1,
    renderer,
    list(),
    list(),
    variables,
    events,
    age,
    parameters,
    timestep
  )

  mockery::expect_args(
    infection_mock,
    1,
    variables,
    events,
    bitten,
    age,
    parameters,
    timestep
  )
})

test_that('simulate_bites integrates eir calculation and mosquito side effects', {
  population <- 4
  timestep <- 5
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(
    list(human_population = population, severe_enabled = TRUE)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  infectivity <- c(.6, 0, .2, .3)
  age <- c(20, 24, 5, 39) * 365

  variables$zeta <- individual::DoubleVariable$new((c(.2, .3, .5, .9)))
  variables$infectivity <- individual::DoubleVariable$new(infectivity)
  variables$mosquito_state <- individual::CategoricalVariable$new(
    c('Sm', 'Pm', 'Im', 'NonExistent'),
    c(rep('Im', 10), rep('Sm', 15), rep('NonExistent', 75))
  )
  variables$species <- individual::CategoricalVariable$new(
    c('All'),
    rep('All', 100)
  )

  lambda_mock <- mockery::mock(.5, .5, .5)
  mosquito_effects_mock <- mockery::mock()
  ode_update <- mockery::mock()

  mockery::stub(simulate_bites, 'effective_biting_rate', lambda_mock)
  mockery::stub(simulate_bites, 'biting_effects_individual', mosquito_effects_mock)
  mockery::stub(simulate_bites, 'rpois', mockery::mock(2))
  mockery::stub(simulate_bites, 'sample.int', mockery::mock(c(2, 3)))
  .pi <- rep(1 / population, population)
  mockery::stub(simulate_bites, 'human_pi', mockery::mock(.pi))
  mockery::stub(simulate_bites, 'mosquito_model_update', ode_update)
  ode_model <- mockery::mock()
  bitten <- simulate_bites(
    renderer,
    list(), #solvers
    list(ode_model), #models
    variables,
    events,
    age,
    parameters,
    timestep
  )

  expect_equal(bitten$to_vector(), c(2, 3))

  f <- parameters$blood_meal_rates[[1]]

  mockery::expect_args(
    lambda_mock,
    1,
    .pi,
    age,
    1,
    list(
      prob_bitten_survives = rep(1, population),
      prob_bitten = rep(1, population),
      prob_repelled = rep(0, population)
    ),
    f,
    1,
    parameters
  )

  effects_args <- mockery::mock_args(mosquito_effects_mock)

  expect_equal(effects_args[[1]][[1]], variables)
  expect_equal(effects_args[[1]][[2]], .55)
  expect_equal(effects_args[[1]][[3]], events)
  expect_equal(effects_args[[1]][[4]], 1)
  expect_equal(effects_args[[1]][[5]]$to_vector(), 11:25)
  expect_equal(effects_args[[1]][[6]]$to_vector(), c(1:10, 11:25))
  expect_equal(effects_args[[1]][[7]], parameters$mum)
  expect_equal(effects_args[[1]][[8]], parameters)
  expect_equal(effects_args[[1]][[9]], timestep)

  mockery::expect_args(ode_update, 1, ode_model, 25, f, parameters$mum)
})

test_that('biting_process integrates mosquito effects and human infection', {
  population <- 4
  timestep <- 5
  parameters <- get_parameters(
    list(human_population = population)
  )

  renderer <- individual::Render$new(5)
  events <- mockery::mock()
  age <- c(20, 24, 5, 39) * 365
  variables <- list(birth = individual::DoubleVariable$new((-age + timestep)))
  lagged_foim <- LaggedValue$new(1, 1)
  lagged_eir <- LaggedValue$new(1, 1)
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)

  biting_process <- create_biting_process(
    renderer,
    solvers,
    models,
    variables,
    events,
    parameters,
    lagged_foim,
    lagged_eir
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
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    lagged_foim,
    lagged_eir,
    1,
    1
  )

  mockery::expect_args(
    infection_mock,
    1,
    variables,
    events,
    bitten,
    age,
    parameters,
    timestep,
    renderer
  )
})

test_that('simulate_bites integrates eir calculation and mosquito side effects', {
  population <- 4
  timestep <- 5
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(
    list(human_population = population)
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

  lambda_mock <- mockery::mock(c(.5, .5, .5, .5))
  mosquito_effects_mock <- mockery::mock()
  eqs_update <- mockery::mock()
  sample_mock <- mockery::mock(c(2, 3))
  pois_mock <- mockery::mock(2)

  mockery::stub(simulate_bites, 'biting_effects_individual', mosquito_effects_mock)
  mockery::stub(simulate_bites, 'rpois', pois_mock)
  mockery::stub(simulate_bites, 'fast_weighted_sample', sample_mock)
  mockery::stub(simulate_bites, 'effective_biting_rates', lambda_mock)
  mockery::stub(simulate_bites, 'aquatic_mosquito_model_update', eqs_update)
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  lagged_foim <- list(LaggedValue$new(12.5, .001))
  lagged_eir <- list(list(LaggedValue$new(12, 10)))
  bitten <- simulate_bites(
    renderer,
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    lagged_foim,
    lagged_eir
  )

  expect_equal(bitten$to_vector(), c(2, 3))

  f <- parameters$blood_meal_rates[[1]]

  effects_args <- mockery::mock_args(mosquito_effects_mock)

  expect_equal(effects_args[[1]][[1]], variables)
  expect_equal(effects_args[[1]][[3]], events)
  expect_equal(effects_args[[1]][[4]], 1)
  expect_equal(effects_args[[1]][[5]]$to_vector(), 11:25)
  expect_equal(effects_args[[1]][[6]]$to_vector(), c(1:10, 11:25))
  expect_equal(effects_args[[1]][[7]], parameters$mum)
  expect_equal(effects_args[[1]][[8]], parameters)
  expect_equal(effects_args[[1]][[9]], timestep)

  mockery::expect_args(eqs_update, 1, models[[1]], 25, f, parameters$mum)
  mockery::expect_args(
    pois_mock,
    1,
    1,
    10 * mean(unique_biting_rate(age, parameters))
  )
  mockery::expect_args(
    sample_mock,
    1,
    2,
    c(.5, .5, .5, .5)
  )
})


test_that('simulate_bites works with mixed populations', {
  population <- 4
  timestep <- 5
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(
    list(human_population = population)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  mock_foim <- mockery::mock(1)
  mock_a <- mockery::mock(.3)

  mockery::stub(simulate_bites, 'calculate_foim', mock_foim)
  mockery::stub(simulate_bites, '.human_blood_meal_rate', mock_a)
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  lagged_foim <- list(LaggedValue$new(12.5, .001), LaggedValue$new(12.5, .01))
  lagged_eir <- list(list(LaggedValue$new(12, 10)), list(LaggedValue$new(12, 10)))
  age <- c(20, 24, 5, 39) * 365
  bitten <- simulate_bites(
    renderer,
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    lagged_foim,
    lagged_eir,
    c(0.2, 0.8),
    2
  )

  mockery::expect_args(mock_foim, 1, .3, c(.001, .01), c(.2, .8))
})

test_that('simulate_bites can halve the mixed transmission for 50% rdt detection', {
  population <- 4
  timestep <- 5
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(
    list(human_population = population, mixing_rdt = TRUE)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)

  mock_rdt <- mockery::mock(.5) # 50%
  mock_foim <- mockery::mock(1)
  mock_pois <- mockery::mock(.3)
  mock_a <- mockery::mock(.3)
  mock_psi <- mockery::mock(rep(1, 4))

  mockery::stub(simulate_bites, 'rdt_detectable', mock_rdt)
  mockery::stub(simulate_bites, 'calculate_foim', mock_foim)
  mockery::stub(simulate_bites, 'rpois', mock_pois)
  mockery::stub(simulate_bites, '.human_blood_meal_rate', mock_a)
  mockery::stub(simulate_bites, 'unique_biting_rate', mock_psi)
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  lagged_foim <- list(LaggedValue$new(12.5, .001), LaggedValue$new(12.5, .01))
  lagged_eir <- list(list(LaggedValue$new(12, 10)), list(LaggedValue$new(12, 10)))
  age <- c(20, 24, 5, 39) * 365
  bitten <- simulate_bites(
    renderer,
    solvers,
    models,
    variables,
    events,
    age,
    parameters,
    timestep,
    lagged_foim,
    lagged_eir,
    c(0.2, 0.8),
    2
  )

  mockery::expect_args(mock_foim, 1, .3, c(.0005, .01), c(.2, .8))
  mockery::expect_args(mock_pois, 1, 1, 9) # = 10 * .2 * .5 + 10 * .8
})

test_that('rdt_detectable adjusts correctly with identity parameters', {
  population <- 4
  parameters <- get_parameters(
    list(human_population = population, rdt_intercept = 0, rdt_coeff = 1)
  )
  variables <- create_variables(parameters)

  mock_detected <- mockery::mock(individual::Bitset$new(4)$insert(c(1, 2))) # 50%

  mockery::stub(rdt_detectable, 'calculate_detected', mock_detected)
  rdt_prev <- rdt_detectable(variables, parameters, 1)

  expect_equal(rdt_prev, 0.5)
})

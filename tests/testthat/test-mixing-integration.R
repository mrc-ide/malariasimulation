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
  mock_a <- mockery::mock(1/3)
  mixed_transmission <- list(
    eir = matrix(c(.2, .8), nrow=1, ncol=2),
    inf = c(.2, .8)
  )
  mock_mixing <- mockery::mock(mixed_transmission, cycle=TRUE)

  mockery::stub(simulate_bites, '.human_blood_meal_rate', mock_a)
  mockery::stub(simulate_bites, 'calculate_foim', mock_foim)
  models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(models, parameters)
  lagged_foim <- LaggedValue$new(12.5, .001)
  lagged_eir <- list(LaggedValue$new(12, 10))
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
    mock_mixing,
    1
  )

  mockery::expect_args(mock_mixing, 1, timestep)
  mockery::expect_args(mock_foim, 1, 1/3, .2)
})

test_that('mixing_fn can return isolated transmission for multiple species', {
  population <- 4
  timestep <- 5
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(
    list(human_population = population)
  )
  parameters <- set_species(
    parameters,
    list(arab_params, fun_params, gamb_params), 
    rep(1/3, 3)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  lagged_infectivity <- list(LaggedValue$new(12.5, .1), LaggedValue$new(12.5, .2))
  lagged_eir <- list(
    list(LaggedValue$new(12, 10), LaggedValue$new(12, 10), LaggedValue$new(12, 10)),
    list(LaggedValue$new(12, 20), LaggedValue$new(12, 20), LaggedValue$new(12, 20))
  )

  mock_rdt <- mockery::mock(.5, cycle = TRUE) # 50%

  mixing_fn <- create_transmission_mixer(
    list(variables, variables),
    list(parameters, parameters),
    lagged_eir,
    lagged_infectivity,
    mixing_tt = 1,
    mixing = list(diag(nrow=2)),
    p_captured_tt = 1,
    p_captured = list(1 - diag(nrow=2)), # full coverage
    p_success = 1
  )

  mockery::stub(mixing_fn, 'rdt_detectable', mock_rdt)

  transmission <- mixing_fn(timestep)
  expect_equal(
    transmission,
    list(
      eir = t(matrix(c(rep(10, 3), rep(20, 3)), nrow=3, ncol=2)),
      inf = c(.1, .2)
    )
  )
})


test_that('mixing_fn can halve the mixed transmission for 50% rdt detection', {
  population <- 4
  timestep <- 5
  renderer <- individual::Render$new(5)
  parameters <- get_parameters(
    list(human_population = population)
  )
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  lagged_infectivity <- list(LaggedValue$new(12.5, .1), LaggedValue$new(12.5, .2))
  lagged_eir <- list(list(LaggedValue$new(12, 10)), list(LaggedValue$new(12, 20)))

  mock_rdt <- mockery::mock(.5, cycle = TRUE) # 50%


  mixing_fn <- create_transmission_mixer(
    list(variables, variables),
    list(parameters, parameters),
    lagged_eir,
    lagged_infectivity,
    mixing_tt = 1,
    mixing = list(matrix(rep(.5, 4), nrow=2, ncol=2)),
    p_captured_tt = 1,
    p_captured = list(1 - diag(nrow=2)), # full coverage
    p_success = 1
  )

  mockery::stub(mixing_fn, 'rdt_detectable', mock_rdt)

  transmission <- mixing_fn(timestep)
  expect_equal(
    transmission,
    list(
      eir = matrix(c(10, 12.5), nrow=2, ncol=1),
      inf = c(.075, .150)
    )
  )
})

test_that('rdt_detectable adjusts correctly with identity parameters', {
  population <- 4
  parameters <- get_parameters(
    list(human_population = population, rdt_intercept = 0, rdt_coeff = 1)
  )
  variables <- list(state = individual::CategoricalVariable$new(
    c('S', 'Tr', 'D', 'A', 'U'),
    c('S', 'Tr', 'A', 'U')
  ))

  rdt_prev <- rdt_detectable(variables, parameters, 1)

  expect_equal(rdt_prev, 0.5)
})

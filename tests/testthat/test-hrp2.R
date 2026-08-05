test_that('set_hrp2_parameters overrides only the supplied values', {
  parameters <- get_parameters()
  updated <- set_hrp2_parameters(parameters, hrp2_shape = 3, hrp2_scale = 40)
  expect_equal(updated$hrp2_shape, 3)
  expect_equal(updated$hrp2_scale, 40)

  unchanged <- set_hrp2_parameters(parameters)
  expect_equal(unchanged$hrp2_shape, parameters$hrp2_shape)
  expect_equal(unchanged$hrp2_scale, parameters$hrp2_scale)
})

test_that('update_hrp2_and_render_infections renders counts and (re)starts the hrp2 clock', {
  timestep <- 10
  population <- 6
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    birth = individual::IntegerVariable$new(rep(0, population)),
    hrp2_infection_time = mock_integer(rep(-1, population))
  )
  renderer <- mock_render(timestep)

  to_D <- individual::Bitset$new(population)$insert(c(1, 2))
  treated <- individual::Bitset$new(population)$insert(c(3))
  to_A <- individual::Bitset$new(population)$insert(c(4))

  update_hrp2_and_render_infections(
    variables,
    renderer,
    parameters,
    timestep,
    to_D,
    treated,
    to_A
  )

  # asymptomatic infections restart the hrp2 clock but are not counted
  # towards the treated/untreated infection totals
  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_untreated_infections',
    2,
    timestep
  )
  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_treated_infections',
    1,
    timestep
  )

  expect_bitset_update(
    variables$hrp2_infection_time$queue_update_mock(),
    timestep,
    c(1, 2, 3, 4)
  )
})

test_that('update_hrp2_and_render_infections does not update hrp2 status when no clinical infections occurred', {
  timestep <- 10
  population <- 6
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    birth = individual::IntegerVariable$new(rep(0, population)),
    hrp2_infection_time = mock_integer(rep(-1, population))
  )
  renderer <- mock_render(timestep)

  to_D <- individual::Bitset$new(population)
  treated <- individual::Bitset$new(population)
  to_A <- individual::Bitset$new(population)

  update_hrp2_and_render_infections(
    variables,
    renderer,
    parameters,
    timestep,
    to_D,
    treated,
    to_A
  )

  mockery::expect_called(variables$hrp2_infection_time$queue_update_mock(), 0)
})

test_that('update_hrp2_and_render_infections restarts the hrp2 clock for new asymptomatic infections without counting them as treated/untreated', {
  timestep <- 10
  population <- 4
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    birth = individual::IntegerVariable$new(rep(0, population)),
    hrp2_infection_time = mock_integer(rep(-1, population))
  )
  renderer <- mock_render(timestep)

  to_D <- individual::Bitset$new(population)
  treated <- individual::Bitset$new(population)
  to_A <- individual::Bitset$new(population)$insert(c(2, 4))

  update_hrp2_and_render_infections(
    variables,
    renderer,
    parameters,
    timestep,
    to_D,
    treated,
    to_A
  )

  mockery::expect_args(renderer$render_mock(), 1, 'n_untreated_infections', 0, timestep)
  mockery::expect_args(renderer$render_mock(), 2, 'n_treated_infections', 0, timestep)

  expect_bitset_update(
    variables$hrp2_infection_time$queue_update_mock(),
    timestep,
    c(2, 4)
  )
})

test_that('update_hrp2_and_render_infections renders age-stratified counts when configured', {
  timestep <- 10
  population <- 4
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_epi_outputs(parameters, hrp2 = c(0, 1825, 100 * 365))

  variables <- list(
    birth = individual::IntegerVariable$new(-c(1, 20, 40, 60) * 365),
    hrp2_infection_time = mock_integer(rep(-1, population))
  )
  renderer <- mock_render(timestep)

  to_D <- individual::Bitset$new(population)$insert(1)
  treated <- individual::Bitset$new(population)$insert(3)
  to_A <- individual::Bitset$new(population)

  update_hrp2_and_render_infections(
    variables,
    renderer,
    parameters,
    timestep,
    to_D,
    treated,
    to_A
  )

  render_names <- vapply(
    mockery::mock_args(renderer$render_mock()),
    function(args) args[[1]],
    character(1)
  )
  expect_in(
    c(
      'n_untreated_infections_0_1824',
      'n_treated_infections_0_1824',
      'n_untreated_infections_1825_36499',
      'n_treated_infections_1825_36499'
    ),
    render_names
  )
})

test_that('create_hrp2_clearance_process clears individuals according to the weibull hazard', {
  timestep <- 30
  population <- 4
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_hrp2_parameters(parameters, hrp2_shape = 2, hrp2_scale = 20)

  # individual 1: never infected (-1, ignored)
  # individual 2: infected at timestep 20 (x = 10)
  # individual 3: infected at timestep 5 (x = 25)
  # individual 4: infected at timestep 29 (x = 1)
  variables <- list(
    hrp2_infection_time = mock_integer(c(-1, 20, 5, 29))
  )

  process <- create_hrp2_clearance_process(variables, parameters)

  cleared_mock <- mockery::mock(c(2))
  mockery::stub(process, 'bernoulli_multi_p', cleared_mock)

  process(timestep)

  positive <- variables$hrp2_infection_time$get_index_of(-1)$not(TRUE)
  expect_equal(positive$to_vector(), c(2, 3, 4))

  x <- timestep - c(20, 5, 29)
  expected_p <- 1 - weibull_survival(x + 1, 2, 20) / weibull_survival(x, 2, 20)
  mockery::expect_args(cleared_mock, 1, expected_p)

  # bernoulli_multi_p was stubbed to return local position 2 of positive (2, 3, 4) -> individual 3
  expect_bitset_update(
    variables$hrp2_infection_time$queue_update_mock(),
    -1,
    c(3)
  )
})

test_that('create_hrp2_clearance_process does nothing when no one is hrp2 positive', {
  timestep <- 30
  population <- 3
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    hrp2_infection_time = mock_integer(rep(-1, population))
  )
  process <- create_hrp2_clearance_process(variables, parameters)
  process(timestep)

  mockery::expect_called(variables$hrp2_infection_time$queue_update_mock(), 0)
})

test_that('create_hrp2_clearance_process guards against numerical survival underflow', {
  timestep <- 1e6
  population <- 1
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_hrp2_parameters(parameters, hrp2_shape = 2, hrp2_scale = 20)

  variables <- list(
    hrp2_infection_time = mock_integer(0)
  )
  process <- create_hrp2_clearance_process(variables, parameters)

  cleared_mock <- mockery::mock(1)
  mockery::stub(process, 'bernoulli_multi_p', cleared_mock)

  process(timestep)

  mockery::expect_args(cleared_mock, 1, 1)
})

test_that('create_hrp2_renderer_process renders the total hrp2 positive count', {
  timestep <- 5
  variables <- list(
    hrp2_infection_time = individual::IntegerVariable$new(c(-1, 1, 2, -1))
  )
  renderer <- mock_render(timestep)
  process <- create_hrp2_renderer_process(renderer, variables)
  process(timestep)

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_hrp2_positive',
    2,
    timestep
  )
})

test_that('create_hrp2_age_renderer_process renders age-stratified hrp2 positive counts', {
  timestep <- 5
  parameters <- get_parameters(list(human_population = 4))
  parameters <- set_epi_outputs(parameters, hrp2 = c(0, 1825, 100 * 365))

  variables <- list(
    birth = individual::IntegerVariable$new(-c(1, 20, 40, 60) * 365),
    hrp2_infection_time = individual::IntegerVariable$new(c(1, -1, 2, 3))
  )
  renderer <- mock_render(timestep)
  process <- create_hrp2_age_renderer_process(variables, parameters, renderer)
  process(timestep)

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_hrp2_positive_0_1824',
    1,
    timestep
  )
  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_hrp2_positive_1825_36499',
    2,
    timestep
  )
})

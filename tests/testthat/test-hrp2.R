test_that('set_hrp2_parameters overrides only the supplied values', {
  parameters <- get_parameters()
  updated <- set_hrp2_parameters(
    parameters,
    hrp2_scale = 40,
    hrp2_asymptomatic_prob = 0.6
  )
  expect_equal(updated$hrp2_scale, 40)
  expect_equal(updated$hrp2_asymptomatic_prob, 0.6)

  unchanged <- set_hrp2_parameters(parameters)
  expect_equal(unchanged$hrp2_scale, parameters$hrp2_scale)
  expect_equal(unchanged$hrp2_asymptomatic_prob, parameters$hrp2_asymptomatic_prob)
})

test_that('update_hrp2_and_render_infections renders counts and marks new infections hrp2 positive', {
  timestep <- 10
  population <- 6
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    birth = individual::IntegerVariable$new(rep(0, population)),
    id = individual::DoubleVariable$new(rep(1, population)),
    hrp2 = mock_integer(rep(0, population))
  )
  renderer <- mock_render(timestep)

  to_D <- individual::Bitset$new(population)$insert(c(1, 2))
  treated <- individual::Bitset$new(population)$insert(c(3))
  to_A <- individual::Bitset$new(population)$insert(c(4))

  # force the asymptomatic infection to be selected as hrp2 positive
  mockery::stub(update_hrp2_and_render_infections, 'bernoulli_multi_p', mockery::mock(1))

  update_hrp2_and_render_infections(
    variables,
    renderer,
    parameters,
    timestep,
    to_D,
    treated,
    to_A
  )

  # asymptomatic infections are marked hrp2 positive but are not counted
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

  # untreated (to_D), treated (Tr) and asymptomatic (to_A) are all marked
  # hrp2 positive in a single update
  expect_bitset_update(
    variables$hrp2$queue_update_mock(),
    1,
    c(1, 2, 3, 4)
  )
})

test_that('update_hrp2_and_render_infections does not update hrp2 status when no clinical infections occurred', {
  timestep <- 10
  population <- 6
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    birth = individual::IntegerVariable$new(rep(0, population)),
    hrp2 = mock_integer(rep(0, population))
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

  mockery::expect_called(variables$hrp2$queue_update_mock(), 0)
})

test_that('update_hrp2_and_render_infections marks new asymptomatic infections hrp2 positive without counting them as treated/untreated', {
  timestep <- 10
  population <- 4
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    birth = individual::IntegerVariable$new(rep(0, population)),
    id = individual::DoubleVariable$new(rep(1, population)),
    hrp2 = mock_integer(rep(0, population))
  )
  renderer <- mock_render(timestep)

  to_D <- individual::Bitset$new(population)
  treated <- individual::Bitset$new(population)
  to_A <- individual::Bitset$new(population)$insert(c(2, 4))

  # force both asymptomatic infections to be selected as hrp2 positive
  mockery::stub(update_hrp2_and_render_infections, 'bernoulli_multi_p', mockery::mock(c(1, 2)))

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
    variables$hrp2$queue_update_mock(),
    1,
    c(2, 4)
  )
})

test_that('hrp2_asymptomatic_probability scales lm detection probability by hrp2_asymptomatic_prob (falciparum)', {
  parameters <- get_parameters()
  parameters <- set_hrp2_parameters(parameters, hrp2_asymptomatic_prob = 0.5)

  detection_mock <- mockery::mock(c(.8, .6, .4))
  mockery::stub(hrp2_asymptomatic_probability, 'probability_of_detection', detection_mock)

  result <- hrp2_asymptomatic_probability(c(100, 200, 300), c(1, 2, 3), parameters)

  mockery::expect_args(detection_mock, 1, c(100, 200, 300), c(1, 2, 3), parameters)
  expect_equal(result, c(.8, .6, .4) * 0.5)
})

test_that('hrp2_asymptomatic_probability assumes vivax asymptomatic infections are always lm-detectable', {
  parameters <- get_parameters(parasite = "vivax")
  parameters <- set_hrp2_parameters(parameters, hrp2_asymptomatic_prob = 0.5)

  result <- hrp2_asymptomatic_probability(c(100, 200), NULL, parameters)

  expect_equal(result, c(0.5, 0.5))
})

test_that('update_hrp2_and_render_infections scales the asymptomatic hrp2 probability by lm detectability (falciparum)', {
  timestep <- 10
  population <- 3
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_hrp2_parameters(parameters, hrp2_asymptomatic_prob = 0.5)

  variables <- list(
    birth = individual::IntegerVariable$new(rep(0, population)),
    id = individual::DoubleVariable$new(c(1, 2, 3)),
    hrp2 = mock_integer(rep(0, population))
  )
  renderer <- mock_render(timestep)

  to_D <- individual::Bitset$new(population)
  treated <- individual::Bitset$new(population)
  to_A <- individual::Bitset$new(population)$insert(c(1, 2, 3))

  probability_mock <- mockery::mock(c(.4, .3, .2))
  mockery::stub(update_hrp2_and_render_infections, 'hrp2_asymptomatic_probability', probability_mock)
  bernoulli_mock <- mockery::mock(integer(0))
  mockery::stub(update_hrp2_and_render_infections, 'bernoulli_multi_p', bernoulli_mock)

  update_hrp2_and_render_infections(
    variables,
    renderer,
    parameters,
    timestep,
    to_D,
    treated,
    to_A
  )

  mockery::expect_args(
    probability_mock,
    1,
    get_age(variables$birth$get_values(to_A), timestep),
    variables$id$get_values(to_A),
    parameters
  )
  mockery::expect_args(bernoulli_mock, 1, c(.4, .3, .2))
})

test_that('update_hrp2_and_render_infections assumes vivax asymptomatic infections are always lm-detectable', {
  timestep <- 10
  population <- 2
  parameters <- get_parameters(parasite = "vivax", overrides = list(human_population = population))
  parameters <- set_hrp2_parameters(parameters, hrp2_asymptomatic_prob = 0.5)

  variables <- list(
    birth = individual::IntegerVariable$new(rep(0, population)),
    hrp2 = mock_integer(rep(0, population))
  )
  renderer <- mock_render(timestep)

  to_D <- individual::Bitset$new(population)
  treated <- individual::Bitset$new(population)
  to_A <- individual::Bitset$new(population)$insert(c(1, 2))

  bernoulli_mock <- mockery::mock(integer(0))
  mockery::stub(update_hrp2_and_render_infections, 'bernoulli_multi_p', bernoulli_mock)

  update_hrp2_and_render_infections(
    variables,
    renderer,
    parameters,
    timestep,
    to_D,
    treated,
    to_A
  )

  mockery::expect_args(bernoulli_mock, 1, c(0.5, 0.5))
})

test_that('update_hrp2_and_render_infections renders age-stratified counts when configured', {
  timestep <- 10
  population <- 4
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_epi_outputs(parameters, hrp2 = c(0, 1825, 100 * 365))

  variables <- list(
    birth = individual::IntegerVariable$new(-c(1, 20, 40, 60) * 365),
    hrp2 = mock_integer(rep(0, population))
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

test_that('create_hrp2_clearance_process only clears hrp2-positive individuals whose infection has cleared (state S)', {
  timestep <- 30
  population <- 5
  parameters <- get_parameters(list(human_population = population))
  parameters <- set_hrp2_parameters(parameters, hrp2_scale = 20)

  # individual 1: hrp2 negative (ignored)
  # individual 2: hrp2 positive, state S (infection cleared - eligible)
  # individual 3: hrp2 positive, state S (infection cleared - eligible)
  # individual 4: hrp2 positive, state S (infection cleared - eligible)
  # individual 5: hrp2 positive, but still infected (state D) - not eligible,
  #   left untouched regardless of the hazard
  variables <- list(
    hrp2 = mock_integer(c(0, 1, 1, 1, 1)),
    state = individual::CategoricalVariable$new(
      c('S', 'A', 'D', 'U', 'Tr'),
      c('A', 'S', 'S', 'S', 'D')
    )
  )

  process <- create_hrp2_clearance_process(variables, parameters)

  # bernoulli_multi_p is called from within process() via bitset_at() - use
  # local_mocked_bindings() so it is intercepted regardless of call depth
  cleared_mock <- mockery::mock(c(2))
  testthat::local_mocked_bindings(bernoulli_multi_p = cleared_mock, .package = "malariasimulation")

  process(timestep)

  # the hazard is only evaluated for the individuals whose infection has
  # cleared (2, 3, 4) - individual 5 (still infected) is excluded
  expected_p <- rate_to_prob(1 / 20)
  mockery::expect_args(cleared_mock, 1, rep(expected_p, 3))

  # bernoulli_multi_p was stubbed to return local position 2 of (2, 3, 4) -> individual 3
  expect_bitset_update(
    variables$hrp2$queue_update_mock(),
    0,
    c(3)
  )
})

test_that('create_hrp2_clearance_process does nothing when no one is hrp2 positive', {
  timestep <- 30
  population <- 3
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    hrp2 = mock_integer(rep(0, population))
  )
  process <- create_hrp2_clearance_process(variables, parameters)
  process(timestep)

  mockery::expect_called(variables$hrp2$queue_update_mock(), 0)
})

test_that('create_hrp2_clearance_process does nothing when hrp2-positive individuals are all still infected', {
  timestep <- 30
  population <- 2
  parameters <- get_parameters(list(human_population = population))

  variables <- list(
    hrp2 = mock_integer(c(1, 1)),
    state = individual::CategoricalVariable$new(
      c('S', 'A', 'D', 'U', 'Tr'),
      c('D', 'A')
    )
  )
  process <- create_hrp2_clearance_process(variables, parameters)
  process(timestep)

  mockery::expect_called(variables$hrp2$queue_update_mock(), 0)
})

test_that('create_hrp2_renderer_process renders the total hrp2 positive count', {
  timestep <- 5
  variables <- list(
    hrp2 = individual::IntegerVariable$new(c(0, 1, 1, 0))
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
    hrp2 = individual::IntegerVariable$new(c(1, 0, 1, 1))
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

test_that('that default rendering works', {
  timestep <- 0
  parameters <- get_parameters()
  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S', 'Tr'),
    c('U', 'A', 'D', 'S')
  )
  birth <- individual::IntegerVariable$new(
    -c(2, 5, 10, 11) * 365
  )
  immunity <- individual::DoubleVariable$new(rep(1, 4))
  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    c('no', 'no', 'yes', 'no')
  )

  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
    is_severe,
    immunity,
    parameters,
    renderer
  )

  mockery::stub(process, 'probability_of_detection', mockery::mock(.5))
  mockery::stub(process, 'bernoulli_multi_p', mockery::mock(1))
  process(timestep)

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_730_3650',
    3,
    timestep
  )
  
  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_detect_730_3650',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'p_detect_730_3650',
    1.5,
    timestep
  )

})

test_that('that default rendering works when no one is in the age range', {
  timestep <- 0
  parameters <- get_parameters()
  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S', 'Tr'),
    rep('S', 4)
  )
  birth <- individual::IntegerVariable$new(
    -c(1, 11, 21, 11) * 365
  )
  immunity <- individual::DoubleVariable$new(rep(1, 4))
  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    rep('no', 4)
  )

  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
    is_severe,
    immunity,
    parameters,
    renderer
  )
  process(timestep)
  mockery::stub(process, 'bernoulli_multi_p', mockery::mock(1))

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_730_3650',
    0,
    timestep
  )
})

test_that('that severe rendering works', {
  timestep <- 0
  year <- 365
  parameters <- get_parameters(list(
    severe_prevalence_rendering_min_ages = c(0, 2) * year,
    severe_prevalence_rendering_max_ages = c(5, 10) * year,
    prevalence_rendering_min_ages = NULL,
    prevalence_rendering_max_ages = NULL
  ))
  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S', 'Tr'),
    c('U', 'D', 'D', 'S')
  )
  birth <- individual::IntegerVariable$new(
    -c(2, 5, 10, 11) * 365
  )
  immunity <- individual::DoubleVariable$new(rep(1, 4))
  is_severe <- individual::CategoricalVariable$new(
    c('yes', 'no'),
    c('no', 'yes', 'no', 'no')
  )
  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
    is_severe,
    immunity,
    parameters,
    renderer
  )
  process(timestep)

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_0_1825',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_severe_0_1825',
    1,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'n_730_3650',
    3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    4,
    'n_severe_730_3650',
    1,
    timestep
  )
})

test_that('that clinical incidence rendering works', {
  timestep <- 0
  year <- 365
  birth <- individual::IntegerVariable$new(
    -c(2, 5, 10, 11) * year
  )

  renderer <- mock_render(1)
  incidence_renderer(
    birth,
    renderer,
    individual::Bitset$new(4)$insert(c(1, 2, 4)),
    individual::Bitset$new(4)$insert(seq(4)),
    c(.1, .2, .3, .4),
    'inc_clinical_',
    c(0, 2) * year,
    c(5, 10) * year,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    1,
    'n_0_1825',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_inc_clinical_0_1825',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'p_inc_clinical_0_1825',
    .3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    4,
    'n_730_3650',
    3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    5,
    'n_inc_clinical_730_3650',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    6,
    'p_inc_clinical_730_3650',
    .6,
    timestep
  )
})

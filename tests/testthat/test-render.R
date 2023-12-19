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

  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
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
    'n_detect_pcr_730_3650',
    3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'n_detect_lm_730_3650',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    4,
    'p_detect_lm_730_3650',
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

  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
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

test_that('that multiple ages are rendered', {
  timestep <- 0

  parameters <- get_parameters(parasite = "vivax")
  parameters <- set_age_outputs(parameters,
                                prevalence = c(2,10,50)*365)

  state <- individual::CategoricalVariable$new(
    c('U', 'A', 'D', 'S', 'Tr'),
    c('U', 'A', 'D', 'S', 'U', 'A', 'D')
  )
  birth <- individual::IntegerVariable$new(
    -c(2, 5, 9, 10, 11, 38, 98) * 365
  )
  immunity <- individual::DoubleVariable$new(rep(1, 4))

  renderer <- mock_render(1)
  process <- create_prevelance_renderer(
    state,
    birth,
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
    'n_730_3649',
    3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    2,
    'n_detect_pcr_730_3649',
    3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    3,
    'n_detect_lm_730_3649',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    4,
    'n_3650_18249',
    3,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    5,
    'n_detect_pcr_3650_18249',
    2,
    timestep
  )

  mockery::expect_args(
    renderer$render_mock(),
    6,
    'n_detect_lm_3650_18249',
    1,
    timestep
  )
})

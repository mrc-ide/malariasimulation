test_that('simulate_infection integrates different types of infection and scheduling', {
  population <- 8
  timestep <- 5
  parameters <- get_parameters(list(
    human_population = population,
    severe_enabled = TRUE
  ))
  events <- create_events(parameters)
  renderer <- mock_render(timestep)

  age <- c(20, 24, 5, 39, 20, 24, 5, 39) * 365
  immunity <- c(.2, .3, .5, .9, .2, .3, .5, .9)
  asymptomatics <- mockery::mock()
  variables <- list(
    ib = individual::DoubleVariable$new(immunity),
    id = individual::DoubleVariable$new(immunity),
    state = list(get_index_of = mockery::mock(asymptomatics))
  )

  total_eir <- 5
  eir <- rep(total_eir / population, population)

  bitten <- individual::Bitset$new(population)$insert(c(1, 3, 5, 7))
  boost_immunity_mock <- mockery::mock()
  infected <- individual::Bitset$new(population)$insert(c(1, 3, 5))
  infection_mock <- mockery::mock(infected)
  clinical <- individual::Bitset$new(population)$insert(c(1, 3))
  clinical_infection_mock <- mockery::mock(clinical)
  severe <- individual::Bitset$new(population)$insert(1)
  severe_infection_mock <- mockery::mock(severe)
  treated <- individual::Bitset$new(population)$insert(3)
  treated_mock <- mockery::mock(treated)
  schedule_mock <- mockery::mock()

  mockery::stub(simulate_infection, 'boost_immunity', boost_immunity_mock)
  mockery::stub(simulate_infection, 'calculate_infections', infection_mock)
  mockery::stub(simulate_infection, 'calculate_clinical_infections', clinical_infection_mock)
  mockery::stub(simulate_infection, 'update_severe_disease', severe_infection_mock)
  mockery::stub(simulate_infection, 'calculate_treated', treated_mock)
  mockery::stub(simulate_infection, 'schedule_infections', schedule_mock)
  simulate_infection(
    variables,
    events,
    bitten,
    age,
    parameters,
    timestep,
    renderer
  )

  mockery::expect_args(
    boost_immunity_mock,
    1,
    variables$ib,
    bitten,
    variables$last_boosted_ib,
    5,
    parameters$ub
  )

  mockery::expect_args(
    infection_mock,
    1,
    variables,
    bitten,
    parameters,
    timestep
  )

  mockery::expect_args(
    clinical_infection_mock,
    1,
    variables,
    infected,
    parameters
  )

  mockery::expect_args(
    severe_infection_mock,
    1,
    timestep,
    clinical,
    variables,
    infected,
    parameters
  )

  mockery::expect_args(
    treated_mock,
    1,
    variables,
    clinical,
    events$recovery,
    events$detection,
    parameters,
    timestep,
    renderer
  )

  mockery::expect_args(
    schedule_mock,
    1,
    events,
    clinical,
    treated,
    infected,
    parameters
  )
})

test_that('calculate_infections works various combinations of drug and vaccination', {
  timestep <- 50
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHA_PQP_params))
  parameters <- set_clinical_treatment(parameters, 2, 1, .5)

  variables <- list(
    state = individual::CategoricalVariable$new(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    drug = individual::DoubleVariable$new(c(1, 2, 0, 0)),
    drug_time = individual::DoubleVariable$new(c(20, 30, -1, -1)),
    rtss_vaccinated = individual::DoubleVariable$new(c(-1, 10, 40, -1)),
    rtss_boosted = individual::DoubleVariable$new(c(-1, 40, -1, -1)),
    rtss_cs = individual::DoubleVariable$new(c(1, .5, 1, 1)),
    rtss_rho = individual::DoubleVariable$new(c(.75, .25, .75, .75)),
    rtss_ds = individual::DoubleVariable$new(c(1, .5, 1, 1)),
    rtss_dl = individual::DoubleVariable$new(c(5, 2, 5, 5)),
    ib = individual::DoubleVariable$new(c(.2, .3, .5, .9))
  )
        
  immunity_mock <- mockery::mock(c(.2, .3, .4))
  weibull_mock <- mockery::mock(.2)
  rtss_antibodies_mock <- mockery::mock(c(2, 3))
  rtss_efficacy_mock <- mockery::mock(c(.2, .3))
  bernoulli_mock <- mockery::mock(2)
  mockery::stub(calculate_infections, 'blood_immunity', immunity_mock)
  mockery::stub(calculate_infections, 'dweibull', weibull_mock)
  mockery::stub(calculate_infections, 'calculate_rtss_antibodies', rtss_antibodies_mock)
  mockery::stub(calculate_infections, 'calculate_rtss_efficacy', rtss_efficacy_mock)
  mockery::stub(calculate_infections, 'bernoulli_multi_p', bernoulli_mock)

  bitten_humans <- individual::Bitset$new(4)$insert(c(1, 2, 3, 4))

  infections <- calculate_infections(
    variables,
    bitten_humans, 
    parameters,
    timestep
  )

  expect_equal(infections$to_vector(), 3)

  mockery::expect_args(immunity_mock, 1, c(.3, .5, .9), parameters)
  mockery::expect_args(
    weibull_mock,
    1,
    20,
    parameters$drug_prophylaxis_shape[[2]],
    parameters$drug_prophylaxis_scale[[2]]
  )
  mockery::expect_args(
    rtss_antibodies_mock,
    1,
    c(10, 10),
    c(.5, 1),
    c(.25, .75),
    c(.5, 1),
    c(2, 5),
    parameters
  )
  mockery::expect_args(
    rtss_efficacy_mock,
    1,
    c(2, 3),
    parameters
  )
  mockery::expect_args(
    bernoulli_mock,
    1,
    c(.2 * .8 * .8, .3 * .7, .4)
  )

})


test_that('calculate_clinical_infections correctly samples clinically infected', {
  timestep <- 5
  parameters <- get_parameters()

  variables <- list(
    ica = individual::DoubleVariable$new(c(.2, .3, .5, .9)),
    icm = individual::DoubleVariable$new(c(.2, .3, .5, .9)),
    id = individual::DoubleVariable$new(c(.2, .3, .5, .9)),
    last_boosted_ica = individual::DoubleVariable$new(c(-1, -1, 1, -1)),
    last_boosted_id = individual::DoubleVariable$new(c(-1, -1, 1, -1))
  )

  immunity_mock <- mockery::mock(c(.2, .3, .4))
  boost_mock <- mockery::mock()
  mockery::stub(calculate_clinical_infections, 'boost_immunity', boost_mock)

  mockery::stub(calculate_clinical_infections, 'clinical_immunity', immunity_mock)
  bernoulli_mock <- mockery::mock(c(1, 3))
  mockery::stub(calculate_clinical_infections, 'bernoulli_multi_p', bernoulli_mock)

  infections <- individual::Bitset$new(4)$insert(c(2, 3, 4))

  clinical_infections <- calculate_clinical_infections(
    variables,
    infections,
    parameters
  )

  expect_equal(clinical_infections$to_vector(), c(2, 4))

  mockery::expect_args(
    immunity_mock,
    1,
    c(.3, .5, .9),
    c(.3, .5, .9),
    parameters
  )

  mockery::expect_args(
    bernoulli_mock,
    1,
    c(.2, .3, .4)
  )
})

test_that('calculate_treated correctly samples treated and updates the drug state', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHA_PQP_params))
  parameters <- set_clinical_treatment(parameters, 1, 1, .25)
  parameters <- set_clinical_treatment(parameters, 2, 1, .25)
  timestep <- 5
  events <- create_events(parameters)
  variables <- list(
    state = list(queue_update = mockery::mock()),
    infectivity = list(queue_update = mockery::mock()),
    drug = list(queue_update = mockery::mock()),
    drug_time = list(queue_update = mockery::mock())
  )

  recovery_mock <- mockery::mock()
  detection_mock <- mockery::mock()
  mockery::stub(calculate_treated, 'recovery$schedule', recovery_mock)
  mockery::stub(calculate_treated, 'detection$schedule', detection_mock)

  seek_treatment <- individual::Bitset$new(4)$insert(c(1, 2, 4))
  mockery::stub(
    calculate_treated,
    'sample_bitset',
    mockery::mock(seek_treatment)
  )
  sample_mock <- mockery::mock(c(2, 1, 1, 1))
  mockery::stub(calculate_treated, 'sample.int', sample_mock)
  bernoulli_mock <- mockery::mock(c(1, 3))
  mockery::stub(calculate_treated, 'bernoulli_multi_p', bernoulli_mock)
  mockery::stub(calculate_treated, 'log_uniform', mockery::mock(c(3, 4)))

  clinical_infections <- individual::Bitset$new(4)
  clinical_infections$insert(c(1, 2, 3, 4))
  calculate_treated(
    variables,
    clinical_infections,
    events$recovery,
    events$detection,
    parameters,
    timestep,
    mock_render(timestep)
  )

  mockery::expect_args(
    sample_mock,
    1,
    2,
    3,
    c(.25, .25),
    TRUE
  )

  mockery::expect_args(
    bernoulli_mock,
    1,
    parameters$drug_efficacy[c(2, 1, 1, 1)]
  )

  expect_bitset_update(variables$state$queue_update, 'Tr', c(1, 4))
  expect_bitset_update(
    variables$infectivity$queue_update,
    parameters$cd * parameters$drug_rel_c[c(2, 1)],
    c(1, 4)
  )
  expect_bitset_update(variables$drug$queue_update, c(2, 1), c(1, 4))
  expect_bitset_update(variables$drug_time$queue_update, 5, c(1, 4))

  expect_bitset_schedule(
    recovery_mock,
    c(1, 4),
    c(3, 4)
  )
  expect_bitset_schedule(
    detection_mock,
    c(1, 4),
    0
  )
})

test_that('schedule_infections correctly schedules new infections', {
  parameters <- get_parameters()

  scheduled <- individual::Bitset$new(20)$insert(c(1, 3, 7, 15))
  clinical_mock <- mockery::mock()
  asym_mock <- mockery::mock()
  detection_mock <- mockery::mock()

  events <- list(
    detection = list(schedule = detection_mock),
    clinical_infection = list(
      schedule = clinical_mock,
      get_scheduled = mockery::mock(individual::Bitset$new(20)$insert(c(1, 3)))
    ),
    asymptomatic_infection = list(
      schedule = asym_mock,
      get_scheduled = mockery::mock(individual::Bitset$new(20)$insert(c(7, 15)))
    ),
    subpatent_infection = list(clear_schedule = mockery::mock()),
    recovery = list(clear_schedule = mockery::mock())
  )

  mockery::stub(
    schedule_infections,
    'log_uniform',
    mockery::mock(
      c(5, 6, 13, 14),
      c(2, 4, 16, 18, 19, 20)
    )
  )

  mockery::stub(
    schedule_infections,
    'bernoulli_multi_p',
    mockery::mock(c(1, 3, 6))
  )

  infections <- individual::Bitset$new(20)$insert(1:20)
  clinical_infections <- individual::Bitset$new(20)$insert(5:15)
  treated <- individual::Bitset$new(20)$insert(7:12)

  schedule_infections(
    events,
    clinical_infections,
    treated,
    infections,
    parameters,
    seq(20),
    rep(1, 20),
    individual::Bitset$new(20)$insert(17)
  )

  expect_bitset_schedule(
    clinical_mock,
    c(5, 6, 13, 14),
    c(5, 6, 13, 14)
  )

  expect_bitset_schedule(
    detection_mock,
    c(5, 6, 13, 14),
    c(5, 6, 13, 14)
  )

  expect_bitset_schedule(
    asym_mock,
    c(2, 4, 16, 18, 19, 20),
    c(2, 4, 16, 18, 19, 20)
  )

  expect_bitset_schedule(
    detection_mock,
    c(2, 4, 16, 18, 19, 20),
    c(2, 4, 16, 18, 19, 20),
    call = 2
  )
})

test_that('prophylaxis is considered for medicated humans', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHA_PQP_params))
  events <- create_events(parameters)
  timestep <- 50

  variables = list(
    state = individual::CategoricalVariable$new(
      c('D', 'S', 'A', 'U', 'Tr'),
      c('D', 'S', 'A', 'U')
    ),
    drug = individual::DoubleVariable$new(c(0, 2, 1, 0)),
    drug_time = individual::DoubleVariable$new(c(-1, 49, 40, -1)),
    rtss_vaccinated = individual::DoubleVariable$new(c(-1, -1, -1, -1)),
    rtss_boosted = individual::DoubleVariable$new(c(-1, -1, -1, -1)),
    ib = individual::DoubleVariable$new(c(.2, .3, .5, .9))
  )

  bitten <- individual::Bitset$new(4)$insert(seq(4))
  m <- mockery::mock(seq(3))
  mockery::stub(calculate_infections, 'bernoulli_multi_p', m)

  calculate_infections(variables, bitten, parameters, timestep)

  expect_equal(
    mockery::mock_args(m)[[1]][[1]],
    c(.590, .384, .590),
    tolerance = 1e-3
  )
})

test_that('boost_immunity respects the delay period', {
  level <- c(2.4, 1.2, 0., 4.)
  immunity <- individual::DoubleVariable$new(level)
  last_boosted <- individual::DoubleVariable$new(c(11, 5, 1, 13))

  level_mock <- mockery::mock()
  mockery::stub(boost_immunity, 'immunity_variable$queue_update', level_mock)

  last_mock <- mockery::mock()
  mockery::stub(boost_immunity, 'last_boosted_variable$queue_update', last_mock)

  index <- individual::Bitset$new(4)$insert(seq(4))
  timestep <- 15
  delay <- 4

  boost_immunity(
    immunity,
    index,
    last_boosted,
    timestep,
    delay
  )

  mockery::expect_args(
    level_mock,
    1,
    c(3.4, 2.2, 1),
    seq(3)
  )

  mockery::expect_args(
    last_mock,
    1,
    15,
    seq(3)
  )
})

test_that('boost_immunity respects the delay period', {
  level <- c(2.4, 1.2, 0., 4., 0.)
  immunity <- individual::DoubleVariable$new(level)
  last_boosted <- individual::DoubleVariable$new(c(11, 5, 1, 13, -1))

  index <- individual::Bitset$new(5)
  index$insert(seq(5))
  timestep <- 15
  delay <- 4

  level_mock <- mockery::mock()
  mockery::stub(boost_immunity, 'immunity_variable$queue_update', level_mock)

  last_mock <- mockery::mock()
  mockery::stub(boost_immunity, 'last_boosted_variable$queue_update', last_mock)

  boost_immunity(
    immunity,
    index,
    last_boosted,
    timestep,
    delay
  )

  mockery::expect_args(
    level_mock,
    1,
    c(3.4, 2.2, 1, 1),
    c(seq(3), 5)
  )

  mockery::expect_args(
    last_mock,
    1,
    15,
    c(seq(3), 5)
  )
})

test_that('boost_immunity does not update when there is no-one to update', {
  level <- c(2.4, 1.2, 0., 4., 0.)
  immunity <- individual::DoubleVariable$new(level)
  last_boosted <- individual::DoubleVariable$new(c(11, 5, 1, 13, -1))

  index <- individual::Bitset$new(5)
  index$insert(seq(5))
  timestep <- 15
  delay <- 4

  level_mock <- mockery::mock()
  mockery::stub(boost_immunity, 'immunity$queue_update', level_mock)

  last_mock <- mockery::mock()
  mockery::stub(boost_immunity, 'last_boosted$queue_update', last_mock)

  boost_immunity(
    immunity,
    index,
    last_boosted,
    timestep,
    delay
  )
  mockery::expect_called(level_mock, 0)
  mockery::expect_called(last_mock, 0)
})




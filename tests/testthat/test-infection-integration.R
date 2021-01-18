test_that('simulate_infection integrates different types of infection and scheduling', {
  population <- 8
  timestep <- 5
  parameters <- get_parameters(list(human_population = population, severe_enabled = TRUE))
  events <- create_events(parameters)

  age <- c(20, 24, 5, 39, 20, 24, 5, 39) * 365
  variables <- list(
    ib = individual::DoubleVariable$new(c(.2, .3, .5, .9, .2, .3, .5, .9))
  )

  total_eir <- 5
  eir <- rep(total_eir / population, population)

  b <- individual::Bitset$new(population)
  b$insert(c(1, 3, 5, 7))
  bernoulli_mock <- mockery::mock(b)
  boost_immunity_mock <- mockery::mock()
  infection_mock <- mockery::mock(c(1, 3, 5))
  clinical_infection_mock <- mockery::mock(c(1, 3))
  severe_infection_mock <- mockery::mock(c(1))
  treated_mock <- mockery::mock(c(3))
  schedule_mock <- mockery::mock()

  mockery::stub(simulate_infection, 'bernoulli_multi_p', bernoulli_mock)
  mockery::stub(simulate_infection, 'boost_immunity', boost_immunity_mock)
  mockery::stub(simulate_infection, 'calculate_infections', infection_mock)
  mockery::stub(simulate_infection, 'calculate_clinical_infections', clinical_infection_mock)
  mockery::stub(simulate_infection, 'update_severe_disease', severe_infection_mock)
  mockery::stub(simulate_infection, 'calculate_treated', treated_mock)
  mockery::stub(simulate_infection, 'schedule_infections', schedule_mock)
  simulate_infection(variables, events, eir, age, parameters, timestep)

  mockery::expect_args(bernoulli_mock, 1, eir)

  mockery::expect_args(
    boost_immunity_mock,
    1,
    variables$ib,
    b,
    c(.2, .5, .2, .5),
    variables$last_boosted_ib,
    5,
    parameters$ub
  )

  mockery::expect_args(
    infection_mock,
    1,
    variables,
    b,
    c(.2, .3, .5, .9, .2, .3, .5, .9),
    parameters,
    timestep
  )

  mockery::expect_args(
    clinical_infection_mock,
    1,
    variables,
    c(1, 3, 5),
    parameters,
    timestep
  )

  mockery::expect_args(
    severe_infection_mock,
    1,
    c(1, 3),
    c(20, 5) * 365,
    variables,
    c(1, 3, 5)
  )

  mockery::expect_args(
    treated_mock,
    1,
    variables,
    c(1, 3),
    events$recovery,
    parameters,
    timestep
  )

  mockery::expect_args(
    schedule_mock,
    1,
    events,
    c(1, 3),
    c(3),
    c(1, 3, 5)
  )
})

test_that('calculate_infections works various combinations of drug and vaccination', {
  timestep <- 50
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, 2, 1)

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
    rtss_dl = individual::DoubleVariable$new(c(5, 2, 5, 5))
  )
        
  immunity_mock <- mockery::mock(c(.2, .3, .4))
  weibull_mock <- mockery::mock(.2)
  rtss_antibodies_mock <- mockery::mock(c(2, 3))
  rtss_efficacy_mock <- mockery::mock(c(.2, .3))
  b <- individual::Bitset$new(4)
  b$insert(2)
  bernoulli_mock <- mockery::mock(b)
  mockery::stub(calculate_infections, 'blood_immunity', immunity_mock)
  mockery::stub(calculate_infections, 'dweibull', weibull_mock)
  mockery::stub(calculate_infections, 'calculate_rtss_antibodies', rtss_antibodies_mock)
  mockery::stub(calculate_infections, 'calculate_rtss_efficacy', rtss_efficacy_mock)
  mockery::stub(calculate_infections, 'bernoulli_multi_p', bernoulli_mock)

  bitten_humans <- individual::Bitset$new(4)
  bitten_humans$insert(c(1, 2, 3, 4))
  infections <- calculate_infections(
    variables,
    bitten_humans, 
    c(.2, .3, .5, .9),
    parameters,
    timestep
  )

  expect_equal(infections$to_vector(), 2)

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
  b <- individual::Bitset$new(4)
  b$insert(c(1, 3))
  bernoulli_mock <- mockery::mock(b)
  mockery::stub(calculate_clinical_infections, 'bernoulli_multi_p', bernoulli_mock)
  infections <- individual::Bitset$new(4)
  infections$insert(c(2, 3, 4))

  clinical_infections <- calculate_clinical_infections(
    variables,
    infections,
    parameters,
    timestep
  )

  expect_equal(clinical_infections$to_vector(), 3)

  mockery::expect_args(
    boost_mock,
    1,
    variables$ica,
    infections,
    c(.3, .5, .9),
    variables$last_boosted_ica,
    5,
    parameters$uc
  )

  mockery::expect_args(
    boost_mock,
    2,
    variables$id,
    infections,
    c(.3, .5, .9),
    variables$last_boosted_id,
    5,
    parameters$ud
  )

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
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, c(1, 2), c(.5, .5))
  timestep <- 5
  events <- create_events(parameters)
  variables <- list(
    state = list(queue_update = mockery::mock()),
    infectivity = list(queue_update = mockery::mock()),
    drug = list(queue_update = mockery::mock()),
    drug_time = list(queue_update = mockery::mock())
  )

  schedule_mock <- mockery::mock()
  mockery::stub(calculate_treated, 'recovery$schedule', schedule_mock)

  seek_treatment <- individual::Bitset$new(4)
  seek_treatment$insert(c(1, 2, 4))
  mockery::stub(
    calculate_treated,
    'sample_bitset',
    mockery::mock(seek_treatment)
  )
  sample_mock <- mockery::mock(c(2, 1, 1, 1))
  mockery::stub(calculate_treated, 'sample.int', sample_mock)
  b <- individual::Bitset$new(4)
  b$insert(c(1, 4))
  bernoulli_mock <- mockery::mock(b)
  mockery::stub(calculate_treated, 'bernoulli_multi_p', bernoulli_mock)
  mockery::stub(calculate_treated, 'log_uniform', mockery::mock(c(3, 4)))

  clinical_infections <- individual::Bitset$new(4)
  clinical_infections$insert(c(1, 2, 3, 4))
  calculate_treated(
    variables,
    clinical_infections,
    events$recovery,
    parameters,
    timestep
  )

  mockery::expect_args(
    sample_mock,
    1,
    2,
    3,
    c(.5, .5),
    TRUE
  )

  mockery::expect_args(
    bernoulli_mock,
    1,
    parameters$drug_efficacy[c(2, 1, 1, 1)]
  )

  expect_equal(
    mockery::mock_args(variables$state$queue_update)[[1]][[1]],
    'Tr'
  )
  expect_equal(
    mockery::mock_args(variables$state$queue_update)[[1]][[2]]$to_vector(),
    c(1, 4)
  )

  expect_equal(
    mockery::mock_args(variables$infectivity$queue_update)[[1]][[1]],
    parameters$cd * parameters$drug_rel_c[c(2, 1)],
  )

  expect_equal(
    mockery::mock_args(variables$infectivity$queue_update)[[1]][[2]]$to_vector(),
    c(1, 4)
  )

  expect_equal(
    mockery::mock_args(variables$drug$queue_update)[[1]][[1]],
    c(2, 1),
  )

  expect_equal(
    mockery::mock_args(variables$drug$queue_update)[[1]][[2]]$to_vector(),
    c(1, 4)
  )

  expect_equal(
    mockery::mock_args(variables$drug_time$queue_update)[[1]][[1]],
    5
  )

  expect_equal(
    mockery::mock_args(variables$drug_time$queue_update)[[1]][[2]]$to_vector(),
    c(1, 4)
  )

  mockery::expect_args(
    schedule_mock,
    1,
    c(1, 4),
    c(3, 4)
  )
})

test_that('schedule_infections correctly schedules new infections', {
  parameters <- get_parameters()

  scheduled <- individual::Bitset$new(20)
  scheduled$insert(c(1, 3, 7, 15))
  clinical_mock <- mockery::mock()
  asym_mock <- mockery::mock()
  all_mock <- mockery::mock(cycle = TRUE)

  events <- list(
    infection = list(
      schedule = all_mock,
      get_scheduled = mockery::mock(scheduled)
    ),
    clinical_infection = list(schedule = clinical_mock),
    asymptomatic_infection = list(schedule = asym_mock),
    subpatent_infection = list(clear_schedule = mockery::mock()),
    recovery = list(clear_schedule = mockery::mock())
  )

  mockery::stub(
    schedule_infections,
    'log_uniform',
    mockery::mock(
      c(5, 6, 13, 14),
      c(2, 4, 16, 17, 18, 19, 20)
    )
  )

  infections <- individual::Bitset$new(20)
  infections$insert(1:20)
  clinical_infections <- individual::Bitset$new(20)
  clinical_infections$insert(5:15)
  treated <- individual::Bitset$new(20)
  treated$insert(7:12)

  schedule_infections(
    events,
    clinical_infections,
    treated,
    infections,
    parameters
  )

  mockery::expect_args(
    clinical_mock,
    1,
    c(5, 6, 13, 14),
    c(5, 6, 13, 14)
  )

  mockery::expect_args(
    all_mock,
    1,
    c(5, 6, 13, 14),
    c(5, 6, 13, 14)
  )

  mockery::expect_args(
    asym_mock,
    1,
    c(2, 4, 16, 17, 18, 19, 20),
    c(2, 4, 16, 17, 18, 19, 20)
  )

  mockery::expect_args(
    all_mock,
    2,
    c(2, 4, 16, 17, 18, 19, 20),
    c(2, 4, 16, 17, 18, 19, 20)
  )
})

test_that('prophylaxis is considered for medicated humans', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  ib = c(.2, .3, .5, .9)
  api <- mock_api(
    list(
      human = list(
        D = c(1),
        S = c(2),
        U = c(3),
        A = c(4),
        drug = c(0, 2, 1, 0),
        drug_time = c(-1, 49, 40, -1)
      )
    ),
    timestep = 50,
    parameters = parameters
  )
  m <- mockery::mock(c(TRUE, TRUE, TRUE, TRUE))

  with_mock(
    'malariasimulation:::bernoulli_multi_p' = m,
    calculate_infections(api, individuals$human, states, variables, seq(4), ib)
  )

  expect_equal(
    mockery::mock_args(m)[[1]][[2]],
    c(0.590, 0.590, 0.384),
    tolerance = 1e-3
  )
})

test_that('boost_immunity respects the delay period', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(
    states,
    variables,
    create_events(),
    parameters
  )

  level <- c(2.4, 1.2, 0., 4.)
  last_boosted <- c(11, 5, 1, 13)
  index <- seq(4)
  timestep <- 15
  delay <- 4

  api <- mock_api(
    list(
      human = list(
        ID = level,
        last_boosted_id = last_boosted
      )
    )
  )

  boost_immunity(
    api,
    individuals$human,
    variables$id,
    index,
    level,
    variables$last_boosted_id,
    timestep,
    delay
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$id,
    c(3.4, 2.2, 1),
    seq(3)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$last_boosted_id,
    15,
    seq(3)
  )
})

test_that('boost_immunity respects the delay period', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(
    states,
    variables,
    create_events(),
    parameters
  )

  level <- c(2.4, 1.2, 0., 4., 0.)
  last_boosted <- c(11, 5, 1, 13, -1)
  index <- seq(5)
  timestep <- 15
  delay <- 4

  api <- mock_api(
    list(
      human = list(
        ID = level,
        last_boosted_id = last_boosted
      )
    )
  )

  boost_immunity(
    api,
    individuals$human,
    variables$id,
    index,
    level,
    variables$last_boosted_id,
    timestep,
    delay
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$id,
    c(3.4, 2.2, 1, 1),
    c(seq(3), 5)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$last_boosted_id,
    15,
    c(seq(3), 5)
  )
})

test_that('boost_immunity does not update when there is no-one to update', {
  parameters <- get_parameters()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(
    states,
    variables,
    create_events(),
    parameters
  )

  level <- c(2.4, 1.2, 0., 4., 0.)
  last_boosted <- c(12, 14, 14, 13, 13)
  index <- seq(5)
  timestep <- 15
  delay <- 4

  api <- mock_api(
    list(
      human = list(
        ID = level,
        last_boosted_id = last_boosted
      )
    )
  )

  boost_immunity(
    api,
    individuals$human,
    variables$id,
    index,
    level,
    variables$last_boosted_id,
    timestep,
    delay
  )

  mockery::expect_called(api$queue_variable_update, 0)
})

test_that('calculate_treated can handle multiple drugs', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, c(1, 2), c(.5, .5))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(human = list(infectivity = c(.1, .2, .3))),
    parameters = parameters,
    timestep = 5
  )

  with_mock(
    'malariasimulation:::bernoulli' = mockery::mock(
      seq(3) # Mock seek treatment
    ),
    'malariasimulation:::bernoulli_multi_p' = mockery::mock(
      c(TRUE, TRUE, FALSE) # Mock drug success
    ),
    'malariasimulation:::log_uniform' = mockery::mock(
      c(10, 12) # Recovery times
    ),
    sample.int = mockery::mock(c(1, 2)),
    calculate_treated(
      api,
      individuals$human,
      states,
      variables,
      seq(3),
      events$recovery
    )
  )

  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$Tr,
    c(1, 2)
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$infectivity,
    parameters$cd * c(AL_params[[2]], DHC_PQP_params[[2]]),
    c(1, 2)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$drug,
    c(1, 2),
    c(1, 2)
  )

  mockery::expect_args(
    api$queue_variable_update,
    3,
    individuals$human,
    variables$drug_time,
    5,
    c(1, 2)
  )

  mockery::expect_args(
    api$schedule,
    1,
    events$recovery,
    c(1, 2),
    c(10, 12)
  )
})

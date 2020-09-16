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
  eir <- rep(
    total_eir / population / length(parameters$blood_meal_rates),
    population
  )
  eir_mock <- mockery::mock(eir, eir, eir)
  mosquito_effects_mock <- mockery::mock()

  mockery::stub(simulate_bites, 'eir', eir_mock)
  mockery::stub(simulate_bites, 'calculate_mosquito_effects', mosquito_effects_mock)
  total_eir <- simulate_bites(api, individuals, states, variables, age, parameters)

  expect_equal(total_eir, eir * 3)

  f <- parameters$blood_meal_rates[[1]]

  mockery::expect_args(
    eir_mock,
    1,
    api,
    individuals$human,
    variables$zeta,
    age,
    1,
    10,
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
    eir,
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

test_that('simulate_infection integrates different types of infection and scheduling', {
  population <- 8
  parameters <- get_parameters(list(human_population = population, severe_enabled = TRUE))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  age <- c(20, 24, 5, 39, 20, 24, 5, 39) * 365

  api <- mock_api(
    list(
      human = list(
        IB = c(.2, .3, .5, .9, .2, .3, .5, .9)
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  total_eir <- 5
  eir <- rep(total_eir / population, population)

  bernoulli_mock <- mockery::mock(c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE))
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
  simulate_infection(api, individuals, states, variables, events, eir, age, parameters)

  mockery::expect_args(bernoulli_mock, 1, 8, eir)

  mockery::expect_args(
    boost_immunity_mock,
    1,
    api,
    individuals$human,
    variables$ib,
    c(1, 3, 5, 7),
    c(.2, .5, .2, .5),
    variables$last_boosted_ib,
    5,
    parameters$ub
  )

  mockery::expect_args(
    infection_mock,
    1,
    api,
    individuals$human,
    states,
    variables,
    c(1, 3, 5, 7),
    c(.2, .3, .5, .9, .2, .3, .5, .9)
  )

  mockery::expect_args(
    clinical_infection_mock,
    1,
    api,
    individuals$human,
    variables,
    c(1, 3, 5)
  )

  mockery::expect_args(
    severe_infection_mock,
    1,
    api,
    c(1, 3),
    c(20, 5) * 365,
    individuals$human,
    variables,
    c(1, 3, 5)
  )

  mockery::expect_args(
    treated_mock,
    1,
    api,
    individuals$human,
    states,
    variables,
    c(1, 3)
  )

  mockery::expect_args(
    schedule_mock,
    1,
    api,
    events,
    c(1, 3),
    c(3),
    c(1, 3, 5)
>>>>>>> 2917837... New feeding cycle:
  )
})

test_that('calculate_infections works various combinations of drug and vaccination', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, 2, 1)
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        D = c(1),
        S = c(2),
        A = c(3),
        U = c(4),
        drug = c(1, 2, 0, 0),
        drug_time = c(20, 30, -1, -1),
        rtss_vaccinated = c(-1, 10, 40, -1),
        rtss_boosted = c(-1, 40, -1, -1),
        rtss_cs = c(1, .5, 1, 1),
        rtss_rho = c(.75, .25, .75, .75),
        rtss_ds = c(1, .5, 1, 1),
        rtss_dl = c(5, 2, 5, 5)
      )
    ),
    timestep = 50,
    parameters = parameters
  )

  immunity_mock <- mockery::mock(c(.2, .3, .4))
  weibull_mock <- mockery::mock(.2)
  rtss_antibodies_mock <- mockery::mock(c(2, 3))
  rtss_efficacy_mock <- mockery::mock(c(.2, .3))
  bernoulli_mock <- mockery::mock(c(FALSE, TRUE, FALSE))
  mockery::stub(calculate_infections, 'blood_immunity', immunity_mock)
  mockery::stub(calculate_infections, 'dweibull', weibull_mock)
  mockery::stub(calculate_infections, 'calculate_rtss_antibodies', rtss_antibodies_mock)
  mockery::stub(calculate_infections, 'calculate_rtss_efficacy', rtss_efficacy_mock)
  mockery::stub(calculate_infections, 'bernoulli_multi_p', bernoulli_mock)

  infections <- calculate_infections(
    api,
    individuals$human,
    states,
    variables,
    c(1, 2, 3, 4), 
    c(.2, .3, .5, .9)
  )

  expect_equal(infections, 3)

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
    3,
    c(.2 * .8 * .8, .3 * .7, .4)
  )

})


test_that('mosquito_effects correctly samples mortalities and infections without interventions', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  infectivity <- c(.6, 0, .2, .3)
  eir <- c(.1, .2, .3, .4)

  api <- mock_api(
    list(),
    parameters = parameters
  )

  f <- parameters$blood_meal_rates[[1]]
  infected <- c(1:25)
  dead <- c(1:15, 51:60, 76:85)
  bernoulli_mock = mockery::mock(infected, dead)
  mockery::stub(calculate_mosquito_effects, 'bernoulli', bernoulli_mock)
  calculate_mosquito_effects(
    api,
    infectivity,
    eir,
    individuals,
    states,
    1,
    1:50,
    1:100,
    1,
    0,
    f,
    parameters
  )

  mockery::expect_args(
    bernoulli_mock,
    1,
    50,
    sum(infectivity * eir)
  )

  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$mosquito,
    states$Pm,
    infected
  )

  mockery::expect_args(
    bernoulli_mock,
    2,
    100,
    parameters$mum
  )

  mockery::expect_args(
    api$queue_state_update,
    2,
    individuals$mosquito,
    states$Unborn,
    dead
  )

})

test_that('calculate_clinical_infections correctly samples clinically infected', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        ICA = c(.2, .3, .5, .9),
        ICM = c(.2, .3, .5, .9),
        ID = c(.2, .3, .5, .9),
        last_boosted_ica = c(-1, -1, 1, -1),
        last_boosted_id = c(-1, -1, 1, -1)
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  immunity_mock <- mockery::mock(c(.2, .3, .4))
  boost_mock <- mockery::mock()
  mockery::stub(calculate_clinical_infections, 'boost_immunity', boost_mock)

  mockery::stub(calculate_clinical_infections, 'clinical_immunity', immunity_mock)
  bernoulli_mock <- mockery::mock(c(TRUE, FALSE, TRUE))
  mockery::stub(calculate_clinical_infections, 'bernoulli_multi_p', bernoulli_mock)
  clinical_infections <- calculate_clinical_infections(
    api,
    individuals$human,
    variables,
    c(2, 3, 4)
  )
  expect_equal(clinical_infections, c(2, 4))

  mockery::expect_args(
    boost_mock,
    1,
    api,
    individuals$human,
    variables$ica,
    c(2, 3, 4),
    c(.3, .5, .9),
    variables$last_boosted_ica,
    5,
    parameters$uc
  )

  mockery::expect_args(
    boost_mock,
    2,
    api,
    individuals$human,
    variables$id,
    c(2, 3, 4),
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
    3,
    c(.2, .3, .4)
  )
})

test_that('calculate_treated correctly samples treated and updates the drug state', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, c(1, 2), c(.5, .5))
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(),
    timestep = 5,
    parameters = parameters
  )

  bernoulli_mock <- mockery::mock(c(1, 2, 4))
  mockery::stub(calculate_treated, 'bernoulli', bernoulli_mock)
  sample_mock <- mockery::mock(c(2, 1, 1))
  mockery::stub(calculate_treated, 'sample.int', sample_mock)
  bernoulli_multi_mock <- mockery::mock(c(TRUE, FALSE, TRUE))
  mockery::stub(calculate_treated, 'bernoulli_multi_p', bernoulli_multi_mock)

  calculate_treated(api, individuals$human, states, variables, c(1, 2, 3, 4))

  mockery::expect_args(
    bernoulli_mock,
    1,
    4,
    parameters$ft
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
    bernoulli_multi_mock,
    1,
    3,
    parameters$drug_efficacy[c(2, 1, 1)]
  )

  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$Tr,
    c(1, 4)
  )
  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$infectivity,
    parameters$cd * parameters$drug_rel_c[c(2, 1)],
    c(1, 4)
  )
  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$drug,
    c(2, 1),
    c(1, 4)
  )
  mockery::expect_args(
    api$queue_variable_update,
    3,
    individuals$human,
    variables$drug_time,
    5,
    c(1, 4)
  )
})

test_that('schedule_infections correctly schedules', {
  skip('not yet')
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(AL_params, DHC_PQP_params))
  parameters <- set_clinical_treatment(parameters, .5, c(1, 2), c(.5, .5))
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

  parameters$severe_enabled = 1

  api <- mock_api(
    list(
      human = list(
        S = c(2),
        U = c(3),
        A = c(4),
        D = c(1),
        birth = -c(20, 24, 5, 39) * 365 + 5,
        IB = c(.2, .3, .5, .9),
        xi = c(.2, .3, .5, .9),
        ICA = c(.2, .3, .5, .9),
        IVA = c(.2, .3, .5, .9),
        ICM = c(.2, .3, .5, .9),
        IVM = c(.2, .3, .5, .9),
        infectivity = c(.6, 0, .2, .3),
        drug = c(0, 0, 0, 0),
        last_boosted_ib = c(-1, -1, 1, -1),
        last_boosted_ica = c(-1, -1, 1, -1),
        last_boosted_iva = c(-1, -1, 1, -1),
        last_boosted_id = c(-1, -1, 1, -1),
        ID = c(.2, .3, .5, .9)
      ),
      mosquito = list(
        Im = 1:100,
        variety = c(rep(1, 25), rep(2, 25), rep(3, 50))
      )
    ),
    timestep = 5,
    parameters = parameters
  )

  bernoulli_mock <- mockery::mock(
    c(TRUE, FALSE, TRUE, FALSE), # bitten
    c(TRUE),                     # infected
    c(TRUE),                     # clinical
    c(FALSE),                    # severe
    c(TRUE)                      # treatment successful
  )

  api$get_scheduled = mockery::mock(4, 1)

  with_mock(
    'malariasimulation:::bernoulli_multi_p' = bernoulli_mock,
    'malariasimulation:::bernoulli' = mockery::mock(1),
    sample.int = mockery::mock(2), # Return DCH drug
    biting_process(api)
  )

  # expect treated individual to go to Tr state
  mockery::expect_args(
    api$queue_state_update,
    1,
    individuals$human,
    states$Tr,
    3
  )

  # expect rel_c to be applied correctly
  mockery::expect_args(
    api$queue_variable_update,
    4,
    individuals$human,
    variables$infectivity,
    parameters$cd * 0.09434,
    3
  )
  
  # expect drug and time to be applied correctly
  mockery::expect_args(
    api$queue_variable_update,
    5,
    individuals$human,
    variables$drug,
    2,
    3
  )
  mockery::expect_args(
    api$queue_variable_update,
    6,
    individuals$human,
    variables$drug_time,
    5,
    3
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
    sample.int = mockery::mock(c(1, 2)),
    calculate_treated(api, individuals$human, states, variables, seq(3))
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
})

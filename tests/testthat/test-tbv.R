test_that('TBV strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_tbv(
    parameters,
    start = 10,
    end = 100,
    frequency = 365,
    ages = c(1, 2, 3, 18),
    coverage = 0.8
  )
  expect_equal(parameters$tbv, TRUE)
  expect_equal(parameters$tbv_start, 10)
  expect_equal(parameters$tbv_end, 100)
  expect_equal(parameters$tbv_ages, c(1, 2, 3, 18))
  expect_equal(parameters$tbv_coverage, .8)
})

test_that('FOIM considers TBA', {
  parameters <- get_parameters()
  parameters <- set_tbv(
    parameters,
    start = 50,
    end = 200,
    frequency = 365,
    ages = c(1, 2, 3, 18),
    coverage = 0.8
  )
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        S = c(1),
        U = c(2),
        A = c(3),
        D = c(4),
        Tr = c(5),
        infectivity = c(0, .1, .15, .5, .3),
        tbv_vaccinated = c(-1, -1, 50, 50, 50)
      )
    ),
    parameters = parameters,
    timestep = 55
  )

  foim_mock = mockery::mock()
  with_mock(
    'malariasimulation:::mosquito_force_of_infection' = foim_mock,
    mosquito_force_of_infection_from_api(
      individuals$human,
      states,
      variables,
      api
    )
  )

  expect_equal(
    mockery::mock_args(foim_mock)[[1]][[4]],
    c(0, .1, 0.275, .352, .325),
    tolerance=1e-3
  )
})

test_that('TBV vaccinations update vaccination time', {
  parameters <- get_parameters()
  parameters <- set_tbv(
    parameters,
    start = 50,
    end = 200,
    frequency = 365,
    ages = c(1, 2, 3, 18),
    coverage = 0.8
  )
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  create_event_based_processes(individuals, states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        birth = -c(8, 2.9, 3.2, 18.4) * 365 + 100,
        rtss_vaccinated = c(-1, -1, 50, 50)
      )
    ),
    parameters = parameters,
    timestep = 100
  )

  with_mock(
    rnorm = mockery::mock(c(0, 0), c(0, 0)),
    'malariasimulation:::bernoulli' = mockery::mock(rep(TRUE, 3)),
    events$tbv_vaccination$listeners[[1]](api, c(2, 3, 4))
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$tbv_vaccinated,
    100,
    c(2, 3, 4)
  )
})

test_that('TBV Antibodies are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_tbv_antibodies(
      c(0, 0, 10, 30),
      parameters
    ),
    c(403.4, 403.4, 116.1, 76.4),
    tolerance = 1e-3
  )
})

test_that('TRAs are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_TRA(c(400, 400, 100, 50), parameters),
    c(0.685, 0.685, 0.466, 0.349),
    tolerance = 1e-3
  )
})


test_that('TBAs are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_TBA(
      c(0.685, 0.685, 0.466, 0.349),
      as.numeric(c(parameters[c('tbv_mt', 'tbv_md', 'tbv_ma', 'tbv_mu')])),
      parameters$tbv_k
    ),
    c(0.685, 0.685, 0.466, 0.349),
    tolerance = 1e-3
  )
})

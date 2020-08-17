test_that('RTS,S strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_rtss(
    parameters,
    start = 10,
    end = 100,
    frequency = 365,
    ages = c(1, 2, 3, 18),
    coverage = 0.8
  )
  expect_equal(parameters$rtss, TRUE)
  expect_equal(parameters$rtss_start, 10)
  expect_equal(parameters$rtss_end, 100)
  expect_equal(parameters$rtss_ages, c(1, 2, 3, 18))
  expect_equal(parameters$rtss_coverage, .8)
})

test_that('Infection considers vaccine efficacy', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)

  api <- mock_api(
    list(
      human = list(
        S = seq(4),
        birth = -c(8, 2.9, 3.2, 18.4) * 365 - 100,
        rtss_vaccinated = c(-1, -1, 50, 50),
        rtss_cs = exp(
          c(rep(parameters$rtss_cs[[1]], 2), rep(parameters$rtss_cs_boost[[1]], 2))
        ),
        rtss_rho = invlogit(
          c(rep(parameters$rtss_rho[[1]], 2), rep(parameters$rtss_rho_boost[[1]], 2))
        ),
        rtss_ds = exp(rep(parameters$rtss_ds[[1]], 4)),
        rtss_dl = exp(rep(parameters$rtss_dl[[1]], 4)),
        drug = c(-1, -1, -1, -1)
      )
    ),
    parameters = parameters,
    timestep = 100
  )

  bernoulli_mock <- mockery::mock(c(TRUE, TRUE, FALSE, FALSE))
  with_mock(
    'malariasimulation:::bernoulli_multi_p' = bernoulli_mock,
    calculate_infections(
      api,
      individuals$human,
      states,
      variables,
      seq(4),
      rep(.2, 4)
    )
  )

  expect_equal(
    mockery::mock_args(bernoulli_mock)[[1]][[2]],
    c(0.590, 0.590, 0.275, 0.275),
    tolerance=1e-3
  )
})

test_that('RTS,S vaccinations update vaccination time and antibody params time', {
  parameters <- get_parameters()
  parameters <- set_rtss(
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
    events$rtss_vaccination$listeners[[1]](api, c(2, 3, 4))
  )

  mockery::expect_args(
    api$queue_variable_update,
    1,
    individuals$human,
    variables$rtss_vaccinated,
    100,
    c(2, 3, 4)
  )

  mockery::expect_args(
    api$queue_variable_update,
    2,
    individuals$human,
    variables$rtss_cs,
    exp(rep(parameters$rtss_cs_boost[[1]], 2)),
    c(3, 4)
  )

  mockery::expect_args(
    api$queue_variable_update,
    3,
    individuals$human,
    variables$rtss_rho,
    invlogit(rep(parameters$rtss_rho_boost[[1]], 2)),
    c(3, 4)
  )
})

test_that('RTS,S Antibodies are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_rtss_antibodies(
      c(0, 0, 10, 30),
      exp(c(6, 6, 5, 5)),
      invlogit(c(2, 2, 1, 1)),
      exp(c(3, 3, 3, 3)),
      exp(c(6, 6, 6, 6)),
      parameters
    ),
    c(403.4, 403.4, 116.1, 76.4),
    tolerance = 1e-3
  )
})

test_that('Efficacies are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_rtss_efficacy(c(400, 400, 100, 50), parameters),
    c(0.685, 0.685, 0.466, 0.349),
    tolerance = 1e-3
  )
})

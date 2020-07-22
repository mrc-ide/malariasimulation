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

  api <- mock_api(list(
    human = list(
      birth = -c(8, 2.9, 3.2, 18.4) * 365 - 100,
      rtss_vaccinated = c(-1, -1, 50, 50),
      rtss_boosted = c(-1, -1, 0, 1),
      drug = c(-1, -1, -1, -1)
    ),
    parameters = parameters,
    timestep = 100
  ))

  bernoulli_mock <- mockery::mock(c(TRUE, TRUE, FALSE, FALSE))
  with_mock(
    'malariasimulation:::bernoulli' = bernoulli_mock,
    calculate_infections(
      api,
      individuals$human,
      states,
      variables,
      seq(4),
      rep(.2, 4)
    )
  )

  mockery::expect_args(
    bernoulli_mock,
    1,
    4,
    c(0, 0, 0, 0)
  )
})

test_that('RTS,S vaccinations update vaccination time and boosted time', {
  parameters <- get_parameters()
  events <- create_events()
  states <- create_states(parameters)
  variables <- create_variables(parameters)
  individuals <- create_individuals(states, variables, events, parameters)
  create_event_based_processes(individuals, states, variables, events, parameters)

  api <- mock_api(list(
    human = list(
      birth = -c(8, 2.9, 3.2, 18.4) * 365 - 100,
      rtss_vaccinated = c(-1, -1, 50, 50),
      rtss_boosted = c(-1, -1, 0, 1)
    ),
    parameters = parameters,
    timestep = 100
  ))

  events$rtss_vaccination$listeners[[1]](api, c(2, 3, 4))

  mockery::expect_args(
    api$queue_variable_update(),
    1,
    individuals$human,
    variables$rtss_vaccinated,
    100,
    c(2, 3, 4)
  )

  mockery::expect_args(
    api$queue_variable_update(),
    2,
    individuals$human,
    variables$rtss_vaccinated,
    c(0, 1, 1),
    c(2, 3, 4)
  )
})

test_that('RTS,S Antibodies are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_rtss_antibodies(
      c(0, 0, 10, 30),
      c(1, 0, 1, 0),
      parameters
    ),
    c(0, 0, 0, 0)
  )
})

test_that('Efficacies are calculated correctly', {
  parameters <- get_parameters()
  expect_equal(
    calculate_rtss_efficacy(c(0, 1, 2, 3, 4), parameters),
    c(0, 0, 0, 0, 0)
  )
})

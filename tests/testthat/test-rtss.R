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
  expect_equal(parameters$rtss_ages, c(1, 2, 3, 18) * 365)
  expect_equal(parameters$rtss_efficacy, .8)
})

test_that('Infection considers vaccine efficacy', {
  parameters <- get_parameters()
  api <- mock_api(list(
    human = list(
      birth = -c(8, 2.9, 3.2, 18.4) * 365 - 100,
      rtss_vaccinated = c(-1, -1, 50, 50),
      rtss_boosted = c(-1, -1, 0, 1)
    ),
    parameters = parameters,
    timestep = 100
  ))

  bernoulli_mock <- mockery::mock(c(TRUE, TRUE, FALSE, FALSE))
  with_mock(
    'malariasimulation:::bernoulli' = bernoulli_mock 
    calculate_infection(api)
  )

  mockery::expect_args(
    bernoulli_mock,
    4,
    c(0, 0, 0, 0)
  )
})

test_that('RTS,S vaccinations update vaccination time and boosted time', {
  api <- mock_api(list(
    human = list(
      birth = -c(8, 2.9, 3.2, 18.4) * 365 - 100,
      rtss_vaccinated = c(-1, -1, 50, 50),
      rtss_boosted = c(-1, -1, 0, 1)
    ),
    parameters = parameters,
    timestep = 100
  ))

  listener(api, c(2, 3, 4))

  mockery::expect_args(
    api$queue_variable_update(),
    individuals$human,
    variables$rtss_vaccinated,
    100,
    c(2, 3, 4)
  )

  mockery::expect_args(
    api$queue_variable_update(),
    individuals$human,
    variables$rtss_vaccinated,
    c(0, 1, 1),
    c(2, 3, 4)
  )
})

test_that('RTS,S Antibodies are calculated correctly on vaccination day', {
  expect_equal(calculate_antibodies(0, parameters), 0)
})

test_that('RTS,S Antibodies are calculated correctly a month later', {
  expect_equal(calculate_antibodies(30, parameters), 0)
})

test_that('RTS,S Antibodies are calculated correctly on boost day', {
  expect_equal(calculate_boosted_antibodies(0, parameters), 0)
})

test_that('RTS,S Antibodies are calculated correctly a month after boosting', {
  expect_equal(calculate_boosted_antibodies(30, parameters), 0)
})

test_that('Impact on infection with no antibodies is 0', {
  expect_equal(rtss_efficacy(0, parameters), 0)
})

test_that('Max antibodies -> Max impact on infection', {
  expect_equal(rtss_efficacy(100, parameters), parameters$vmax)
})

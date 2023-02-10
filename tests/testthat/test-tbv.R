test_that('TBV strategy parameterisation works', {
  parameters <- get_parameters()
  parameters <- set_tbv(
    parameters,
    timesteps = 10,
    coverages = 0.8,
    ages = c(1, 2, 3, 18)
  )
  expect_equal(parameters$tbv, TRUE)
  expect_equal(parameters$tbv_timesteps, 10)
  expect_equal(parameters$tbv_coverages, .8)
  expect_equal(parameters$tbv_ages, c(1, 2, 3, 18))
  
  expect_error(
    parameters <- set_tbv(
      parameters,
      timesteps = 10,
      coverages = -1,
      ages = c(1, 2, 3, 18)
    ), "all(coverages >= 0) & all(coverages <= 1) is not TRUE",
    fixed = TRUE
  )
  
  expect_error(
    parameters <- set_tbv(
      parameters,
      timesteps = 10,
      coverages = 1.5,
      ages = c(1, 2, 3, 18)
    ), "all(coverages >= 0) & all(coverages <= 1) is not TRUE",
    fixed = TRUE
  )
})

test_that('TBV scheduler works on first timestep', {
  parameters <- get_parameters()
  parameters <- set_tbv(
    parameters,
    timesteps = 10,
    coverages = 0.8,
    ages = c(1, 2, 3, 18)
  )
  output <- run_simulation(10, parameters)
  expect_equal(output$n_vaccinated_tbv[1:9], as.numeric(rep(NA, 9)))
  expect_gt(output$n_vaccinated_tbv[[10]], 0)
})

test_that('TBV antibodies are calculated correctly', {
  tau <- 22
  rho <- .7
  ds <- 45
  dl <- 591
  t <- c( 0, 0, 10, 30 )
  expected <- c( 22, 22, 19.7, 16.1 )
  expect_equal(calculate_tbv_antibodies(t, tau, rho, ds, dl), expected, tolerance=1e-1)
})

test_that('TRA is calculated correctly', {
	mu <- 12.63
	gamma1 <- 2.5
	gamma2 <- .06
	antibodies <- c( 22, 22, 19.7, 5 );
	expected <- c( 0.985, 0.985, 0.981, 0.622 );
	expect_equal(calculate_TRA(mu, gamma1, gamma2, antibodies), expected, tolerance=1e-3)
})

test_that('TBA is calculated correctly', {
	tra <- c( 0.685, 0.685, 0.466, 0.349 )
	mx <- c( 35, 46.7, 3.6, .8 )
	k <- .9
	expected <- c( 0.06379336, 0.04997884, 0.16019879, 0.22684080 )
	expect_equal(calculate_TBA(mx, k, tra), expected, tolerance=1e-3)
})

test_that('account_for_tbv integrates tbv functions correctly', {
  timestep <- 55
  parameters <- get_parameters(list(human_population = 5))
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  variables$state <- individual::CategoricalVariable$new(
    c('S', 'A', 'D', 'U', 'Tr'),
    c('S', 'U', 'A', 'D', 'Tr')
  )
  variables$tbv_vaccinated <- individual::DoubleVariable$new(
    c(-1, -1, 50, 50, 50)
  )
  infectivity <- c( 0, .1, .15, .5, .3 )
  expected <- c( 0.0, 0.1, 0.0112133312, 0.2226002769, 0.1137443524 )
  expect_equal(
    account_for_tbv(
      timestep,
      infectivity,
      variables,
      parameters
    ),
    expected
  )
})

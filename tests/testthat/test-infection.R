test_that('probability_bitten returns correct values', {
  age <- c(0, 5, 30)
  xi <-  c(1.8, 2., .5)
  infectious_variants <- c(rep(1, 3), rep(3, 2))

  days_per_timestep <- 1
  parameters <- list(
    rho   = .85,
    a0    = 8 * 365 / days_per_timestep,
    av1   = .92,
    av2   = .74,
    av3   = .94
  )
  expect_equal(
    probability_bitten(
      age,
      xi,
      infectious_variants,
      parameters
    ),
    c(.4739, .5964, .0409),
    tolerance=1e-4
  )
})

test_that('blood_immunity returns correct values', {
  ib <-  c(0., 1., 10.)

  parameters <- list(
    b0    = 0.590076,
    b1    = 0.5,
    ib0   = 43.8787,
    kb    = 2.15506
  )
  expect_equal(
    blood_immunity(
      ib,
      parameters
    ),
    c(.590, .590, .578),
    tolerance=1e-3
  )
})

test_that('clinical immunity returns correct values', {
  acquired_immunity <-  c(0., 1., 10.)
  maternal_immunity <-  c(.97, .5, .3)
  parameters <- list(
    phi0  = .0749886,
    phi1  = .0001191,
    ic0     = 18.02366,
    kc      = 2.36949
  )
  expect_equal(
    clinical_immunity(acquired_immunity, maternal_immunity, parameters),
    c(.0749, .0748, .0593),
    tolerance=1e-4
  )
})

test_that('severe immunity returns correct values', {
  acquired_immunity <-  c(0., 1., 6.)
  maternal_immunity <-  c(.97, .5, .3)
  age <- c(0, 5, 30)
  parameters <- list(
    theta0  = .0749886,
    theta1  = .0001191,
    kv      = 2.00048,
    fv0     = 0.141195,
    av      = 2493.41,
    gammav  = 2.91282,
    iv0     = 1.09629
  )
  expect_equal(
    severe_immunity(age, acquired_immunity, maternal_immunity, parameters),
    c(0.0675, 0.0593, 0.0132),
    tolerance=1e-4
  )
})

test_that('create_infectivity_frame constructs the correct dataframe', {
  age <- c(0, 5, 30)
  xi <-  c(1.8, 2., .5)
  subset_to_param <-  list(list(c(1, 2), .5), list(c(3), .2))
  expect_mapequal(
    create_infectivity_frame(age, xi, subset_to_param),
    data.frame(age = c(0, 5, 30), xi = c(1.8, 2., .5), infectivity=c(.5, .5, .2))
  )
})

test_that('create_infectivity_frame survives empty subsets', {
  age <- c(0, 5, 30)
  xi <-  c(1.8, 2., .5)
  subset_to_param <-  list(list(c(1,2), .5), list(c(), .2))
  expect_mapequal(
    create_infectivity_frame(age, xi, subset_to_param),
    data.frame(age = c(0, 5), xi = c(1.8, 2.), infectivity=c(.5, .5))
  )
})

test_that('mosquito_force_of_infection returns correct values', {
  days_per_timestep <- 1
  parameters <- list(
    av1   = .92,
    av2   = .74,
    av3   = .94,
    rho   = .85,
    a0    = 8 * 365 / days_per_timestep
  )
  human_frame <- data.frame(
    age = c(0, 5, 30),
    xi = c(1.8, 2., .5),
    infectivity=c(.5, .5, .2)
  )
  v <- c(1, 1, 1, 2, 3, 3)
  expect_equal(
    mosquito_force_of_infection(v, human_frame, parameters),
    c(.426, .426, .426, .343, .436, .436),
    tolerance=1e-3
  )
})

test_that('remove_scheduled removes indecies of already scheduled infections', {
  subset <- c(1, 2, 4)
  current_schedule <- c(-1, 4, 6, 9, -1, 3)
  timestep <- 5
  expect_equal(
    remove_scheduled(subset, timestep, current_schedule),
    c(1, 2)
  )
})

test_that('boost_acquired_immunity respects the delay period', {
  level <- c(2.4, 1.2, 0., 4.)
  last_boosted <- c(11, 5, 1, 13)
  timestep <- 15
  delay <- 4
  expect_equal(
    boost_acquired_immunity(level, last_boosted, timestep, delay),
    c(2.4, 2.2, 1., 4.)
  )
})

test_that('boost_acquired_immunity works with when never boosted', {
  level <- c(2.4, 1.2, 0., 4.)
  last_boosted <- c(2, 5, 1, -1)
  timestep <- 6
  delay <- 10
  expect_equal(
    boost_acquired_immunity(level, last_boosted, timestep, delay),
    c(2.4, 1.2, 0., 5.)
  )
})

test_that('asymptomatic_infectivity returns the correct values', {
  age <- c(0, 5, 30, 6)
  immunity <- c(2.4, 1.2, 0., 4.)
  days_per_timestep <- 1
  parameters <- list(
    cu    = 0.00062,
    cd    = 0.068,
    gamma1= 1.82425,
    fd0   = 0.007055,
    ad    = 21.9 * 365 * days_per_timestep,
    gammad= 4.8183,
    d1    = 1,
    dmin  = 0.161, #NOTE: what should this be?
    id0   = 1.577533,
    kd    = .476614
  )
  expect_equal(
    asymptomatic_infectivity(age, immunity, parameters),
    c(.207, .208, .208, .207),
    tolerance=1e-3
  )
})

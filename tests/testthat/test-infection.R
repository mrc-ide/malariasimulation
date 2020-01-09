test_that('force_of_infection returns correct values', {
  age <- c(0, 5, 30)
  xi <-  c(1.8, 2., .5)
  ib <-  c(0., 1., 6.)
  infectious_variants <- c(rep(1, 3), rep(3, 2))

  timestep_to_day <- 1
  parameters <- list(
    rho   = .85,
    a0    = 8 * 365 * timestep_to_day,
    b0    = 0.590076,
    b1    = 0.5,
    ib0   = 43.8787,
    kb    = 2.15506,
    av1   = .92,
    av2   = .74,
    av3   = .94
  )
  expect_equal(
    force_of_infection(
      age,
      xi,
      infectious_variants,
      ib,
      parameters
    ),
    c(.2369, .2982, .0205),
    tolerance=1e-4
  )
})

test_that('immunity returns correct values', {
  acquired_immunity <-  c(0., 1., 6.)
  maternal_immunity <-  c(.97, .5, .3)
  parameters <- list(
    theta0  = .0749886,
    theta1  = .0001191,
    ic0     = 18.02366,
    kc      = 2.36949
  )
  expect_equal(
    immunity(acquired_immunity, maternal_immunity, parameters),
    c(.0751, .0752, .0812),
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
  timestep_to_day <- 1
  parameters <- list(
    av1   = .92,
    av2   = .74,
    av3   = .94,
    rho   = .85,
    a0    = 8 * 365 * timestep_to_day
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

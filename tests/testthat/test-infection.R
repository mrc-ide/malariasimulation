test_that('vector infectivity in the IBM returns correct values', {
  infectious_variants <- c(rep(1, 3), rep(3, 2))
  parameters <- list(
    blood_meal_rates = c(.92, .74, .94)
  )
  expect_equal(
    vector_infectivity_ibm(
      infectious_variants,
      parameters
    ),
    4.64
  )
})

test_that('vector infectivity in the ODE returns correct values', {
  infectious_variants <- c(rep(1, 3), rep(3, 2))
  ode_states <- mockery::mock(
    c(0, 0, 0, 0, 0, 3),
    c(0, 0, 0, 0, 0, 0),
    c(0, 0, 0, 0, 0, 2)
  )
  parameters <- list(
    blood_meal_rates = c(.92, .74, .94)
  )
  mockery::stub(vector_infectivity_ode, 'mosquito_model_get_states', ode_states)
  expect_equal(
    vector_infectivity_ode(
      list(mockery::mock(), mockery::mock(), mockery::mock()),
      parameters
    ),
    4.64
  )
})

test_that('eir returns correct values', {
  age <- c(0, 5, 30) * 365
  xi <-  c(1.8, 2., .5)

  days_per_timestep <- 1
  parameters <- list(
    rho   = .85,
    a0    = 8 * 365 / days_per_timestep
  )
  infectivity <- 4.64

  expect_equal(
    eir(
      age,
      xi,
      c(1, 3),
      infectivity,
      parameters
    ),
    c(.677, 1.229),
    tolerance=1e-3
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

test_that('mosquito_force_of_infection returns correct values', {
  days_per_timestep <- 1
  parameters <- list(
    blood_meal_rates = c(.92, .74, .94),
    rho   = .85,
    a0    = 8 * 365 / days_per_timestep
  )
  age <- c(0, 5, 30) * 365
  xi <- c(1.8, 2., .5)
  infectivity <- c(.5, .5, .2)
  v <- c(1, 1, 1, 2, 3, 3)
  expect_equal(
    mosquito_force_of_infection(
      v,
      age,
      xi,
      infectivity,
      c(1, 2, 3),
      parameters
    ),
    c(.387, .387, .387, .311, .395, .395),
    tolerance=1e-3
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
    ad    = 21.9 * 365 / days_per_timestep,
    gammad= 4.8183,
    d1    = .2,
    id0   = 1.577533,
    kd    = .476614
  )
  expect_equal(
    asymptomatic_infectivity(age, immunity, parameters),
    c(0.06716227, 0.06732774, 0.06800000, 0.06664776),
    tolerance=1e-3
  )
})

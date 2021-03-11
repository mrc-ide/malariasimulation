test_that('total_M and EIR functions are consistent with equilibrium EIR', {
  EIR <- 5
  jamie_parameters <- malariaEquilibrium::load_parameter_set()
  h_eq <- malariaEquilibrium::human_equilibrium_no_het(
    EIR,
    0,
    jamie_parameters,
    0:99
  )
  foim <- sum(h_eq[,'inf'] * h_eq[,'psi'])
  population <- 100000
  parameters <- get_parameters(c(
    translate_jamie(remove_unused_jamie(jamie_parameters)),
    list(
      init_foim = foim,
      species = 'all',
      species_proportions = 1,
      human_population = population
    )
  ))
  parameters <- parameterise_mosquito_equilibrium(parameters, EIR)
  m_eq <- initial_mosquito_counts(parameters, foim)

  #set up arguments for EIR calculation
  variables <- create_variables(parameters)
  age <- get_age(variables$birth$get_values(), 0)
  zeta <- rep(1, length(age))
  variables$zeta <- individual::DoubleVariable$new(zeta)
  p_bitten <- prob_bitten(0, variables, 1, parameters)
  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(zeta, psi)
  Z <- average_p_repelled(p_bitten$prob_repelled, .pi, parameters$Q0)
  W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, parameters$Q0)
  f <- blood_meal_rate(1, Z, parameters)
  estimated_eir <- m_eq[[6]] * effective_biting_rate(
    .pi,
    age,
    1,
    p_bitten,
    f,
    W,
    parameters
  )
  expect_equal(
    mean(estimated_eir) * 365,
    EIR,
    tolerance = 1
  )
})

test_that('total_M and EIR functions are consistent with equilibrium EIR (with het)', {
  population <- 100000

  EIR <- 50

  jamie_parameters <- malariaEquilibrium::load_parameter_set()

  h_eq <- malariaEquilibrium::human_equilibrium(
    EIR,
    0,
    jamie_parameters,
    0:99
  )
  foim <- h_eq$FOIM

  parameters <- get_parameters(c(
    translate_jamie(remove_unused_jamie(jamie_parameters)),
    list(
      init_foim = foim,
      species = 'all',
      species_proportions = 1,
      human_population = population
    )
  ))
  parameters <- parameterise_mosquito_equilibrium(parameters, EIR)
  m_eq <- initial_mosquito_counts(parameters, foim)

  #set up arguments for EIR calculation
  variables <- create_variables(parameters)
  age <- get_age(variables$birth$get_values(), 0)
  zeta <- variables$zeta$get_values()
  p_bitten <- prob_bitten(0, variables, 1, parameters)
  psi <- unique_biting_rate(age, parameters)
  .pi <- human_pi(zeta, psi)
  Z <- average_p_repelled(p_bitten$prob_repelled, .pi, parameters$Q0)
  W <- average_p_successful(p_bitten$prob_bitten_survives, .pi, parameters$Q0)
  f <- blood_meal_rate(1, Z, parameters)
  estimated_eir <- m_eq[[6]] * effective_biting_rate(
    .pi,
    age,
    1,
    p_bitten,
    f,
    W,
    parameters
  )

  expect_equal(
    mean(estimated_eir) * 365,
    EIR,
    tolerance = 1
  )
})

test_that('mosquito_limit is set to 1 for 0 EIR', {
  parameters <- parameterise_mosquito_equilibrium(get_parameters(), 0)
  expect_equal(parameters$mosquito_limit, 1)
})

test_that('mosquito_limit is set to a sensible level', {
  EIRs <- c(5, 50, 1000)

  seasonalities <- list(
    list(
      g0 = 7.2,
      g = c(-10, 4.06, -1.08),
      h = c(-4.6, 3.59, -1.32)
    ),
    list(
      g0 = 7.2,
      g = c(-10, 4.06, -1.08),
      h = c(-4.6, 3.59, -1.32)
    )
  )

  for (EIR in EIRs) {
    for (seasonality in seasonalities) {
      parameters <- get_parameters(c(
        seasonality,
        model_seasonality = TRUE,
        init_foim = .1
      ))
      parameters <- parameterise_mosquito_equilibrium(parameters, EIR)
      run_simulation(365, parameters)
    }
  }
  expect_true(TRUE)
})

test_that('mosquito_effects correctly samples mortalities and infections without interventions', {
  parameters <- get_parameters()
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  infectivity <- c(.6, 0, .2, .3)
  lambda <- c(.1, .2, .3, .4)
  f <- parameters$blood_meal_rates[[1]]
  infected <- individual::Bitset$new(100)$insert(1:25)
  dead <- individual::Bitset$new(100)$insert(c(1:15, 51:60, 76:85))
  bernoulli_mock = mockery::mock(infected, dead)
  mockery::stub(calculate_mosquito_effects, 'sample_bitset', bernoulli_mock)
  events$mosquito_infection <- mock_event(events$mosquito_infection)
  events$mosquito_death <- mock_event(events$mosquito_death)
  variables$mosquito_state <- mock_category(
    c('Sm', 'Pm', 'Im', 'NonExistent'),
    rep('Sm', 100)
  )
  calculate_mosquito_effects(
    variables,
    infectivity,
    lambda,
    events,
    1,
    individual::Bitset$new(100)$insert(1:50),
    individual::Bitset$new(100)$insert(1:100),
    1,
    0,
    f,
    parameters,
    individual::Render$new(1),
    timestep = 1
  )

  mockery::expect_args(
    variables$mosquito_state$queue_update,
    1,
    'Pm',
    infected
  )

  expect_equal(
    mockery::mock_args(bernoulli_mock)[[2]][[2]],
    parameters$mum
  )

  mockery::expect_args(
    events$mosquito_death$schedule,
    1,
    dead,
    0
  )
})

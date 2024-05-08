test_that('total_M and EIR functions are consistent with equilibrium EIR', {
  skip_on_ci()
  expected_EIR <- c(1, 5, 10, 100)
  population <- 100000
  parameters <- get_parameters(
    list(
      human_population = population,
      enable_heterogeneity = FALSE
    )
  )
  actual_EIR <- vnapply(expected_EIR, function(EIR) {
    parameters <- set_equilibrium(parameters, EIR)

    #set up arguments for EIR calculation
    variables <- create_variables(parameters)
    vector_models <- parameterise_mosquito_models(parameters, 1)
    solvers <- parameterise_solvers(vector_models, parameters)
    estimated_eir <- calculate_eir(1, solvers, variables, parameters, 0)
    age <- get_age(variables$birth$get_values(), 0)
    psi <- unique_biting_rate(age, parameters)
    omega <- mean(psi)
    mean(estimated_eir) / omega / population * 365
  })

  expect_equal(
    actual_EIR * 0.766,
    # 0.766 is the population average psi under an exponential pop distribution
    # There's possibly a better way of calculating this, depending on what exactly
    # we're testing here, e.g., we could remove lines 19-21 and omega from line 22
    expected_EIR,
    tolerance = 1e-2
  )
})

test_that('total_M and EIR functions are consistent with equilibrium EIR (with het)', {
  skip_on_ci()
  population <- 100000

  expected_EIR <- c(1, 5, 10, 100)

  parameters <- get_parameters(list(human_population = population))

  actual_EIR <- vnapply(expected_EIR, function(EIR) {
    parameters <- set_equilibrium(parameters, EIR)

    variables <- create_variables(parameters)
    vector_models <- parameterise_mosquito_models(parameters, 1)
    solvers <- parameterise_solvers(vector_models, parameters)
    estimated_eir <- calculate_eir(1, solvers, variables, parameters, 0)
    age <- get_age(variables$birth$get_values(), 0)
    psi <- unique_biting_rate(age, parameters)
    omega <- mean(psi)
    mean(estimated_eir) / omega / population * 365
  })

  expect_equal(
    actual_EIR * 0.766,
    # 0.766 is the population average psi under an exponential pop distribution
    # There's possibly a better way of calculating this, depending on what exactly
    # we're testing here, e.g., we could remove lines 49-51 and omega from line 52
    expected_EIR,
    tolerance = 1e-2
  )
})

test_that('FOIM is consistent with equilibrium', {
  skip_on_ci()
  population <- 100000

  EIRs <- c(1, 5, 10, 100)

  eq_params <- malariaEquilibrium::load_parameter_set()

  ages <- 0:999 / 10
  expected_foim <- vnapply(
    EIRs,
    function(EIR) {
      malariaEquilibrium::human_equilibrium(
        EIR,
        0,
        eq_params,
        ages
      )$FOIM
    }
  )

  actual_foim <- vnapply(
    EIRs,
    function(EIR) {
      parameters <- get_parameters(c(list(human_population = population)))
      parameters <- set_equilibrium(parameters, EIR)

      variables <- create_variables(parameters)
      vector_models <- parameterise_mosquito_models(parameters, 1)
      solvers <- parameterise_solvers(vector_models, parameters)
      a <- human_blood_meal_rate(1, variables, parameters, 0)
      age <- get_age(variables$birth$get_values(), 0)
      psi <- unique_biting_rate(age, parameters)
      zeta <- variables$zeta$get_values()
      .pi <- human_pi(psi, zeta)
      calculate_foim(a, sum(.pi * variables$infectivity$get_values()), 1.)
    }
  )
  expect_equal(
    expected_foim,
    actual_foim,
    tolerance = 1e-3
  )
})

test_that('phi is consistent with equilibrium at high EIR (no het)', {
  skip_on_ci()
  population <- 1e5

  EIR <- 100
  ages <- 0:999 / 10

  parameters <- get_parameters(list(
    human_population = population,
    enable_heterogeneity = FALSE
  ))
  parameters <- set_equilibrium(parameters, EIR)
  eq_params <- malariaEquilibrium::load_parameter_set()

  eq <-  malariaEquilibrium::human_equilibrium_no_het(
    EIR,
    0,
    eq_params,
    ages
  )

  variables <- create_variables(parameters)
  expect_equal(
    mean(
      clinical_immunity(
        variables$ica$get_values(),
        variables$icm$get_values(),
        parameters
      )
    ),
    weighted.mean(eq[,'phi'], eq[,'prop']),
    tolerance=1e-2
  )
})

test_that('phi is consistent with equilibrium at high EIR', {
  skip_on_ci()
  population <- 1e5

  EIR <- 100
  ages <- 0:999 / 10

  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, EIR)
  eq_params <- malariaEquilibrium::load_parameter_set()

  eq <-  malariaEquilibrium::human_equilibrium(
    EIR,
    0,
    eq_params,
    ages
  )

  variables <- create_variables(parameters)
  expect_equal(
    mean(
      clinical_immunity(
        variables$ica$get_values(),
        variables$icm$get_values(),
        parameters
      )
    ),
    weighted.mean(eq$states[,'phi'], eq$states[,'prop']),
    tolerance=1e-2
  )
})

test_that('mosquito_limit is set to 0 for 0 EIR', {
  parameters <- parameterise_mosquito_equilibrium(get_parameters(list(
    individual_mosquitoes = TRUE
  )), 0)
  expect_equal(parameters$mosquito_limit, 0)
})

test_that('mosquito_limit is set to a sensible level', {
  EIRs <- 5

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
      run_simulation(5, parameters)
    }
  }
  expect_true(TRUE)
})

test_that('mosquito_effects correctly samples mortalities and infections without interventions', {
  parameters <- get_parameters()
  events <- create_events(parameters)
  variables <- create_variables(parameters)
  infected <- individual::Bitset$new(100)$insert(1:25)
  dead <- individual::Bitset$new(100)$insert(c(1:15, 51:60, 76:85))
  bernoulli_mock = mockery::mock(infected, dead)
  mockery::stub(biting_effects_individual, 'sample_bitset', bernoulli_mock)
  events$mosquito_infection <- mock_event(events$mosquito_infection)
  events$mosquito_death <- mock_event(events$mosquito_death)
  variables$mosquito_state <- mock_category(
    c('Sm', 'Pm', 'Im', 'NonExistent'),
    rep('Sm', 100)
  )
  biting_effects_individual(
    variables,
    .5,
    events,
    1,
    individual::Bitset$new(100)$insert(1:50),
    individual::Bitset$new(100)$insert(1:100),
    parameters$mum,
    parameters,
    timestep = 1
  )

  mockery::expect_args(
    variables$mosquito_state$queue_update_mock(),
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

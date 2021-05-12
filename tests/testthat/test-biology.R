test_that('total_M and EIR functions are consistent with equilibrium EIR', {
  EIR <- 5
  population <- 100000
  parameters <- get_parameters(
    list(
      human_population = population,
      enable_heterogeneity = FALSE
    )
  )
  parameters <- set_equilibrium(parameters, EIR)

  #set up arguments for EIR calculation
  variables <- create_variables(parameters)
  vector_models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(vector_models, parameters)
  omega <- sum(unique_biting_rate(-variables$birth$get_values(), parameters))
  estimated_eir <- calculate_eir(1, solvers, variables, parameters, 0) / omega
  expect_equal(
    mean(estimated_eir) * 365,
    EIR,
    tolerance = 1e-2
  )
})

test_that('total_M and EIR functions are consistent with equilibrium EIR (with het)', {
  population <- 100000

  EIR <- 50

  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, EIR)

  variables <- create_variables(parameters)
  vector_models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(vector_models, parameters)
  omega <- sum(unique_biting_rate(-variables$birth$get_values(), parameters))
  estimated_eir <- calculate_eir(1, solvers, variables, parameters, 0) / omega
  expect_equal(
    mean(estimated_eir) * 365,
    EIR,
    tolerance = 1
  )
})

test_that('FOIM is consistent with equilibrium', {
  population <- 1e5

  EIR <- 5

  eq_params <- malariaEquilibrium::load_parameter_set()

  ages <- 0:999 / 10
  foim <- malariaEquilibrium::human_equilibrium(
    EIR,
    0,
    eq_params,
    ages
  )$FOIM

  parameters <- get_parameters(c(list(human_population = population)))
  parameters <- set_equilibrium(parameters, EIR)

  variables <- create_variables(parameters)
  vector_models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(vector_models, parameters)
  lambda <- effective_biting_rates(1, variables, parameters, 0)
  expect_equal(
    calculate_foim(variables$infectivity$get_values(), lambda),
    foim,
    tolerance=1e-2
  )
})

test_that('FOI is consistent with equilibrium', {
  population <- 1e5

  EIR <- 5

  eq_params <- malariaEquilibrium::load_parameter_set()

  ages <- 0:999 / 10
  eq <- malariaEquilibrium::human_equilibrium(
    EIR,
    0,
    eq_params,
    ages
  )
  foi <- weighted.mean(eq$states[,'FOI'], eq$states[,'prop'])

  parameters <- get_parameters(list(human_population = population))
  parameters <- set_equilibrium(parameters, EIR)

  variables <- create_variables(parameters)
  lambda <- effective_biting_rates(1, variables, parameters, 0)
  ib <- blood_immunity(variables$ib$get_values(), parameters)
  expect_equal(
    mean(lambda * ib),
    foi,
    tolerance=1e-2
  )
})

test_that('phi is consistent with equilibrium at high EIR (no het)', {
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

test_that('mosquito_limit is set to 1 for 0 EIR', {
  parameters <- parameterise_mosquito_equilibrium(get_parameters(), 0)
  expect_equal(parameters$mosquito_limit, 1)
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

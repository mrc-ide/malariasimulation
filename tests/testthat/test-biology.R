test_that('total_M and EIR functions are consistent with equilibrium EIR', {
  EIR <- 5
  eq_params <- malariaEquilibrium::load_parameter_set()
  omega <- 1 - eq_params$rho * eq_params$eta / (eq_params$eta + 1/eq_params$a0)
  alpha <- eq_params$f * eq_params$Q0
  states <- malariaEquilibrium::human_equilibrium_no_het(
    EIR,
    0,
    eq_params,
    0:99
  )
  eq <- list(
    states = states,
    FOIM = sum(states[,'inf'] * states[,'psi']) * alpha / omega
  )
  population <- 100000
  parameters <- get_parameters(c(
    translate_jamie(remove_unused_jamie(eq_params)),
    list(
      init_foim = eq$FOIM,
      species = 'all',
      species_proportions = 1,
      human_population = population,
      enable_heterogeneity = FALSE
    )
  ))
  parameters <- parameterise_mosquito_equilibrium(parameters, EIR)
  m_eq <- initial_mosquito_counts(parameters, 1, eq$FOIM, parameters$total_M)

  #set up arguments for EIR calculation
  variables <- create_variables(parameters)
  vector_models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(vector_models, parameters)
  estimated_eir <- calculate_eir(1, solvers, variables, parameters, 0)
  expect_equal(
    mean(estimated_eir) * 365,
    EIR
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
  m_eq <- initial_mosquito_counts(parameters, 1, foim, parameters$total_M)

  #set up arguments for EIR calculation
  variables <- create_variables(parameters)
  vector_models <- parameterise_mosquito_models(parameters)
  solvers <- parameterise_solvers(vector_models, parameters)
  estimated_eir <- calculate_eir(1, solvers, variables, parameters, 0)
  expect_equal(
    mean(estimated_eir) * 365,
    EIR
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

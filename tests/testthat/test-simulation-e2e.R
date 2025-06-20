test_that('Simulation runs for a few timesteps', {
  sim <- run_simulation(100)
  expect_equal(nrow(sim), 100)
})

test_that('run_metapop_simulation fails with incorrect mixing matrix', {
  population <- 4
  timestep <- 5
  parameters <- get_parameters(list(human_population = population))
  paramlist <- list(parameters, parameters)
  # incorrect params
  mixing <- matrix(c(1, 1), nrow = 1, ncol = 2)
  p_captured <- diag(nrow=2)
  expect_error(
    run_metapop_simulation(
      timesteps,
      parameters,
      NULL,
      mixing_tt = 1,
      export_mixing = list(mixing),
      import_mixing = list(mixing),
      p_captured_tt = 1,
      p_captured = list(diag(nrow=2)),
      p_success = 1
    )
  )
})

test_that('run_metapop_simulation integrates two models correctly', {
  population <- 4
  timesteps <- 5
  parameters <- get_parameters(list(human_population = population))
  parametersets <- list(parameters, parameters)
  mixing <- diag(nrow = 2)
  p_captured <- 1 - diag(nrow = 2)
  
  outputs <- run_metapop_simulation(
    timesteps,
    parametersets,
    NULL,
    mixing_tt = 1,
    export_mixing = list(mixing),
    import_mixing = list(mixing),
    p_captured_tt = 1,
    p_captured = list(p_captured),
    p_success = 1
  )
  expect_equal(length(outputs), 2)
  expect_equal(nrow(outputs[[1]]), 5)
  expect_equal(nrow(outputs[[2]]), 5)
})

test_that("run_simulation_with_repetitions() runs successfully without correlations specfied", {
  
  # Load the default parameters:
  parameters <- get_parameters()
  
  # Specify a number of repetitions to run:
  reps <- 3
  
  # Check the run_simulation_with_repetitions function runs without any correlation object specified:
  testthat::expect_no_error(
    simulation <- run_simulation_with_repetitions(
      timesteps = 10,
      parameters = parameters, 
      parallel = F,
      repetitions = reps))
  
  # Check that the repetitions present in the output matches expectations:
  expect_identical(object = unique(simulation$repetition), expected = 1:reps)
  
})

test_that("run_simulation_with_repetitions() runs successfully with correlations specfied", {
  
  # Load the default model parameters:
  parameters <- get_parameters()
  
  # Set some vaccination strategy
  parameters <- set_mass_pev(
    parameters,
    profile = rtss_profile,
    timesteps = 1,
    coverages = .9,
    min_wait = 0,
    min_ages = 100,
    max_ages = 1000,
    booster_spacing = numeric(0),
    booster_coverage = numeric(0),
    booster_profile = NULL
  )
  
  # Set some smc strategy
  parameters <- set_drugs(parameters, list(SP_AQ_params))
  parameters <- set_smc(
    parameters,
    drug = 1,
    timesteps = 100,
    coverages = .9,
    min_age = 100,
    max_age = 1000
  )
  
  # Correlate the vaccination and smc targets
  correlations <- get_correlation_parameters(parameters)
  correlations$inter_intervention_rho('pev', 'smc', 1)
  
  # Correlate the rounds of smc
  correlations$inter_round_rho('smc', 1)
  
  # Specify a number of repetitions to simulate:
  reps <- 2
  
  # Run the simulation without any correlation object specified:
  testthat::expect_no_error(
    output <- run_simulation_with_repetitions(
      timesteps = 10,
      parameters = parameters, 
      correlations = correlations,
      parallel = F,
      repetitions = reps
    )
  )
  
  # Check that the repetitions present in the output matches expectations:
  expect_identical(object = unique(output$repetition), expected =  1:reps)
  
})

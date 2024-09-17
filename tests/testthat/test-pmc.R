test_that("pmc parameterisation works", {
  p <- get_parameters()
  
  expect_false(p$pmc)
  expect_equal(p$pmc_ages, NULL)
  expect_equal(p$pmc_coverages, NULL)
  expect_equal(p$pmc_timesteps, NULL)
  
  
  p <- set_drugs(
    parameters = p,
    drugs = list(SP_AQ_params_falciparum))
  p <- set_pmc(
    parameters = p,
    drug = 1,
    timesteps = c(5, 10),
    coverages = c(0.5, 1),
    ages = c(30, 60, 90)
  )
  
  expect_error(
    parameters <- set_pmc(
      parameters = p,
      drug = 1,
      timesteps = c(5, 10),
      coverages = c(-0.5, 0.5),
      ages = c(30, 60, 90)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE)
  
  expect_error(
    parameters <- set_pmc(
      parameters = p,
      drug = 1,
      timesteps = c(5, 10),
      coverages = c(0.5, 1.5),
      ages = c(30, 60, 90)
    ), "all(coverages >= 0) && all(coverages <= 1) is not TRUE",
    fixed = TRUE)
  
  expect_true(p$pmc)
  expect_equal(p$pmc_ages, c(30, 60, 90))
  expect_equal(p$pmc_coverages, c(0.5, 1))
  expect_equal(p$pmc_timesteps, c(5, 10))
  expect_equal(p$pmc_drug, 1)
})

test_that("pmc gives drugs to correct ages", {

  p <- get_parameters(list(human_population = 6))
  p <- set_drugs(
    parameters = p,
    drugs = list(SP_AQ_params_falciparum))

  p <- set_pmc(
    parameters = p,
    drug = 1,
    timesteps = 10,
    coverages = 1,
    ages = c(30, 60, 90)
  )
  timestep <- 10
  events <- create_events(p)
  renderer <- individual:::Render$new(timestep)
  variables <- create_variables(p)
  variables$birth <- individual::IntegerVariable$new(
    -c(10, 30, 60, 89, 90, 36500) + 10
  )
  variables$state <- mock_category(
    c('D', 'S', 'A', 'U', 'Tr'),
    c('D', 'S', 'A', 'U', 'D', 'S')
  )
  variables$drug <- mock_integer(rep(0, 6))
  variables$drug_time <- mock_integer(rep(-1, 6))
  mockery::stub(sample_intervention, 'bernoulli', mockery::mock(c(TRUE, TRUE, TRUE)))
  local_mocked_bindings(bernoulli_multi_p = mockery::mock(1:3))
  local_mocked_bindings(calculate_asymptomatic_detectable = mockery::mock(individual::Bitset$new(6)$insert(3)))
  
  process <- create_pmc_process(
    variables = variables,
    events = events, 
    parameters = p,
    renderer = renderer,
    correlations = get_correlation_parameters(p),
    coverages = p$pmc_coverages,
    timesteps = p$pmc_timesteps,
    drug = p$pmc_drug
  )

  process(timestep)
  
  # Three treatments given
  expect_equal(renderer$to_dataframe(),
               data.frame(timestep = 1:10,
                          n_pmc_treated = c(rep(0, 9), 3),
                          n_pmc_drug_efficacy_failures = c(rep(0, 10)),
                          n_pmc_successfully_treated = c(rep(0, 9), 3)))
  
  # Individuals 3 and 5, are correct age and in D or A states
  expect_bitset_update(
    variables$state$queue_update_mock(),
    'Tr',
    c(3, 5),
    1
  )
  # Individual 2 is correct age and in S state
  expect_bitset_update(
    variables$state$queue_update_mock(),
    'S',
    2,
    2
  )
  # Drug is recorded as given
  expect_bitset_update(
    variables$drug$queue_update_mock(),
    1,
    c(2, 3, 5)
  )
  # Drug time is recorded
  expect_bitset_update(
    variables$drug_time$queue_update_mock(),
    10,
    c(2, 3, 5)
  )
})

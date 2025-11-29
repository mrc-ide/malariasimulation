test_that("set_nmf sets rates correctly", {
  p <- get_parameters()
  p <- set_nmf(p, ages = c(365, 10 * 365), rates = c(0.1, 0.05))
  expect_equal(p$nmf_ages, c(365, 10 * 365))
  expect_equal(p$nmf_rates, c(0.1, 0.05))
})

test_that("set_nmf validates inputs", {
  p <- get_parameters()
  expect_error(set_nmf(p, ages = c(1, 1), rates = c(0.5)), "length")
  expect_error(set_nmf(p, ages = c(10, 5), rates = c(0.1, 0.1)), "diff")
  expect_error(set_nmf(p, ages = 1, rates = -0.1), ">= 0")
})

test_that("non malarial fevers + patent infection = individuals treated", {
  p <- get_parameters(list(human_population = 5))
  p <- set_drugs(p, list(AL_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = 0, coverages = 1)
  p <- set_nmf(p, ages = c(200 * 365), rates = c(Inf))

  vars <- create_variables(p)
  states <- c('S', 'D', 'A', 'U', 'Tr')
  vars$state <- individual::CategoricalVariable$new(states, c('S', 'D', 'A', 'U', 'Tr'))
  renderer <- individual::Render$new(1)
  nmf <- individual::Bitset$new(p$human_population)

  local_mocked_bindings(
    rdt_detectable = mockery::mock(rep(1, 5))
  )

  nmf_process <- create_nmf_process(vars, p, renderer, nmf)
  nmf_process(1)

  falciparum_infection_outcome_process(
    1,
    individual::Bitset$new(p$human_population),
    nmf,
    vars,
    renderer,
    p
  )

  df <- renderer$to_dataframe()
  # All individuals get a NMF
  expect_equal(df$n_nmf[[1]], 5)
  # Only those in patent states get ("D", "A") detected
  expect_equal(df$n_nmf_malaria_detected [[1]], 2)
})

test_that("non malarial fevers do not occur when set to 0", {
  p <- get_parameters(list(human_population = 5))
  p <- set_drugs(p, list(AL_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = 0, coverages = 1)
  p <- set_nmf(p, ages = c(200 * 365), rates = c(0))
  
  vars <- create_variables(p)
  states <- c('S', 'D', 'A', 'U', 'Tr')
  vars$state <- individual::CategoricalVariable$new(states, c('S', 'D', 'A', 'U', 'Tr'))
  renderer <- individual::Render$new(1)
  nmf <- individual::Bitset$new(p$human_population)
  
  local_mocked_bindings(
    rdt_detectable = mockery::mock(rep(1, 5))
  )
  
  nmf_process <- create_nmf_process(vars, p, renderer, nmf)
  nmf_process(1)
  
  falciparum_infection_outcome_process(
    1,
    individual::Bitset$new(p$human_population),
    nmf,
    vars,
    renderer,
    p
  )
  
  df <- renderer$to_dataframe()
  # No individuals get a NMF
  expect_equal(df$n_nmf[[1]], 0)
})

test_that("non malarial fevers are overriden if occuring alongside a malaria fever", {
  p <- get_parameters(list(human_population = 5))
  p <- set_drugs(p, list(AL_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = 0, coverages = 1)
  p <- set_nmf(p, ages = c(200 * 365), rates = c(Inf))
  
  vars <- create_variables(p)
  states <- c('S', 'D', 'A', 'U', 'Tr')
  vars$state <- individual::CategoricalVariable$new(states, c('A', 'A', 'A', 'A', 'A'))
  renderer <- individual::Render$new(1)
  nmf <- individual::Bitset$new(p$human_population)
  
  local_mocked_bindings(
    rdt_detectable = mockery::mock(rep(1, 5))
  )
  
  nmf_process <- create_nmf_process(vars, p, renderer, nmf)
  nmf_process(1)
  infected_humans <- nmf$copy()
  
  local_mocked_bindings(
    rdt_detectable = mockery::mock(rep(1, 5)),
    calculate_clinical_infections = mockery::mock(infected_humans)
  )
  
  falciparum_infection_outcome_process(
    1,
    infected_humans,
    nmf,
    vars,
    renderer,
    p
  )
  
  df <- renderer$to_dataframe()
  # All individuals get a NMF
  expect_equal(df$n_nmf[[1]], 5)
  expect_equal(df$n_infections[[1]], 5)
  expect_equal(df$n_treated [[1]], 5)
  # NMFs are overriden by clinical malaria infections
  expect_equal(df$n_nmf_malaria_detected [[1]], 0)
})

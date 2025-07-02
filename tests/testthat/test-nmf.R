test_that('set_nmf sets rates correctly', {
  p <- get_parameters()
  p <- set_nmf(p, ages = c(365, 10*365), rates = c(0.1, 0.05))
  expect_equal(p$nmf_ages, c(365, 10*365))
  expect_equal(p$nmf_rates, c(0.1, 0.05))
})

test_that('set_nmf validates inputs', {
  p <- get_parameters()
  expect_error(set_nmf(p, ages = c(1,1), rates = c(0.5)), 'length')
  expect_error(set_nmf(p, ages = c(10,5), rates = c(0.1,0.1)), 'diff')
  expect_error(set_nmf(p, ages = 1, rates = -0.1), 'all\(rates >= 0\)')
})


test_that('non malarial fevers treat individuals', {
  p <- get_parameters(list(human_population = 3))
  p <- set_drugs(p, list(SP_AQ_params))
  p <- set_clinical_treatment(p, drug = 1, timesteps = 0, coverages = 1)
  p <- set_nmf(p, ages = c(1000), rates = c(1))

  vars <- create_variables(p)
  renderer <- individual::Render$new(1)
  nmf <- individual::Bitset$new(p$human_population)

  local_mocked_bindings(
    rdt_detectable = function(variables, parameters, timestep) 1
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
  expect_equal(df$n_nmf[[1]], 3)
  expect_equal(df$n_treated_nmf[[1]], 3)
  expect_equal(vars$state$get_values(), rep('Tr', 3))
  expect_equal(vars$nmf_count$get_values(), rep(1L,3))
})

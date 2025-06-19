test_that('default params can be translated to equilibrium', {
  malariaEquilibrium::human_equilibrium(
    5,
    0,
    translate_parameters(get_parameters()),
    0:100
  )

  expect_true(TRUE)
})

test_that('set_equilibrium does not duplicate eq_params when called multiple times on same parameter list', {
  
  parameters <- get_parameters()
  parameters |>
    set_equilibrium(
      init_EIR = 16
    ) -> parameters
  expect_equal(object = sum(names(parameters) == "eq_params"), expected = 1)
  
  parameters$dt <- 0.001
  parameters |>
    set_equilibrium(init_EIR = 16) -> parameters
  expect_equal(object = sum(names(parameters) == "eq_params"), expected = 1)
  
})

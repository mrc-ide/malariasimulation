test_that('default params can be translated to equilibrium', {
  malariaEquilibrium::human_equilibrium(
    5,
    0,
    translate_parameters(get_parameters()),
    0:100
  )

  expect_true(TRUE)
})


test_that('You can set time varying clinical treatment', {
  parameters <- get_parameters()
  parameters <- set_drugs(parameters, list(DHC_PQP_params, AL_params))
  parameters <- set_clinical_treatment(
    parameters,
    drug = 1,
    timesteps = c(50, 100),
    coverages = c(.2, .5)
  )
  parameters <- set_clinical_treatment(
    parameters,
    drug = 2,
    timesteps = c(25, 75),
    coverages = c(.3, .4)
  )

  expect_equal(
    get_treatment_coverages(parameters, 10),
    c(0, 0)
  )

  expect_equal(
    get_treatment_coverages(parameters, 25),
    c(0, .3)
  )

  expect_equal(
    get_treatment_coverages(parameters, 30),
    c(0, .3)
  )

  expect_equal(
    get_treatment_coverages(parameters, 50),
    c(.2, .3)
  )

  expect_equal(
    get_treatment_coverages(parameters, 80),
    c(.2, .4)
  )
})

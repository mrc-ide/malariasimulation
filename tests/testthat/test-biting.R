test_that('compartmental always gives positive infectious', {
  solver_states <- rep(0, length(ODE_INDICES) + length(ADULT_ODE_INDICES))
  solver_states[[ADULT_ODE_INDICES['Im']]] <- -1e-10
  expect_gte(calculate_infectious_compartmental(solver_states), 0)
})

test_that('gonotrophic_cycle cannot be negative', {
  params <- get_parameters()
  vparams <- gamb_params
  vparams$blood_meal_rates <- 5
  expect_error(set_species(params, list(vparams), 1))
})

test_that('calculate_foim mixes correctly', {
  expect_equal( # check 1 case
    calculate_foim(a=0.5, infectivity=6, mixing=1),
    3
  )
  expect_equal( # check 2 cases
    calculate_foim(a=0.5, infectivity=c(6, 4), mixing=c(.5, .5)),
    2.5
  )
})

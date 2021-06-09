test_that('create_variables allows empty species', {
  params <- malariasimulation::get_parameters()
  params <- malariasimulation::set_species(
    params,
    species=list(
      malariasimulation::gamb_params,
      malariasimulation::arab_params,
      malariasimulation::fun_params
    ),
    proportions=c(1,0,0)
  )
  variables <- create_variables(params)
  expect_equal(variables$species$get_size_of('gamb'), params$mosquito_limit)
})

test_that('create_variables allows multiple species', {
  params <- malariasimulation::get_parameters()
  params <- malariasimulation::set_species(
    params,
    species=list(
      malariasimulation::gamb_params,
      malariasimulation::arab_params
    ),
    proportions=c(.9, .1)
  )
  variables <- create_variables(params)
  expect_equal(
    variables$species$get_size_of('arab'),
    params$total_M * .1
  )
})

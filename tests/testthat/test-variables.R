test_that('create_variables allows empty species', {
  params <- get_parameters()
  params <- set_species(
    params,
    species=list(
      gamb_params,
      arab_params,
      fun_params
    ),
    proportions=c(1,0,0)
  )
  variables <- create_variables(params)
  expect_equal(variables$species$get_size_of('gamb'), params$mosquito_limit)
})

test_that('create_variables allows multiple species', {
  params <- get_parameters()
  params <- set_species(
    params,
    species=list(
      gamb_params,
      arab_params
    ),
    proportions=c(.9, .1)
  )
  variables <- create_variables(params)
  expect_equal(
    variables$species$get_size_of('arab'),
    params$total_M * .1
  )
})

test_that('create_variables allows multiple species w different total_M', {
  params <- get_parameters()
  params <- set_species(
    params,
    species=list(
      gamb_params,
      arab_params
    ),
    proportions=c(.9, .1)
  )
  params <- parameterise_total_M(params, 1000)
  variables <- create_variables(params)
  expect_equal(
    variables$species$get_size_of('arab'),
    params$total_M * .1
  )
})

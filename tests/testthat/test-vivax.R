test_that('Test model runs using falciparum', {
  parameters <- get_parameters()
  sim_res <- run_simulation(timesteps = 100, parameters = parameters)
  expect_equal(nrow(sim_res), 100)
})

test_that('Test vivax model runs', {
  parameters <- get_parameters()
  vivax_parameters <- set_parasite(parameters = parameters, parasite = "vivax")
  sim_res <- run_simulation(timesteps = 100, parameters = vivax_parameters)
  expect_equal(nrow(sim_res), 100)
})

## Test subpatent progression functions
test_that('Test subpatent duration function works ', {
  parameters <- get_parameters(overrides = list(s_proportion = 0,
                                                d_proportion = 0,
                                                a_proportion = 0,
                                                u_proportion = 1,
                                                t_proportion = 0))
  vivax_parameters <- set_parasite(parameters = parameters, parasite = "vivax")
  variables <- create_variables(vivax_parameters)
  index <- variables$state$get_index_of("U")
  expect_equal(object = subpatent_duration(parameters = vivax_parameters, variables = variables, index = index), 
               expected = rep(vivax_parameters$du_max,100))
  
  ## Change initial values of ID, and IDM, check they are the same
  variables$id <- individual::DoubleVariable$new(1:100)
  ID_durations <- subpatent_duration(parameters = vivax_parameters, variables = variables, index = index)
  
  variables$id <- individual::DoubleVariable$new(rep(0,100))
  variables$idm <- individual::DoubleVariable$new(1:100)
  IDM_durations <- subpatent_duration(parameters = vivax_parameters, variables = variables, index = index)
  
  expect_equal(object = ID_durations, expected = IDM_durations)
  
  ## Check convergence to min_du at high immunity
  variables$id <- individual::DoubleVariable$new(rep(100000,100))
  variables$idm <- individual::DoubleVariable$new(rep(0,100))
  expect_equal(object = subpatent_duration(parameters = vivax_parameters, variables = variables, index = index),
               expected = rep(vivax_parameters$du_min,100))
})

test_that('Test idm rendered creates output', {
  parameters <- get_parameters()
  vivax_parameters <- set_parasite(parameters = parameters, parasite = "vivax")
  sim_res <- run_simulation(timesteps = 200, parameters = vivax_parameters)
  expect_vector(object = sim_res$idm_mean, size = 200) 
})


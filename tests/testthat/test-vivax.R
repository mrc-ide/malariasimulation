test_that('Test falciparum switch produces same', {
  parameters_def <- get_parameters()
  parameters_fal <- get_parameters(parasite = "falciparum")
  expect_identical(parameters_def, parameters_fal)
})

test_that('Test vivax model runs', {
  vivax_parameters <- get_parameters(parasite = "vivax")
  sim_res <- run_simulation(timesteps = 100, parameters = vivax_parameters)
  expect_equal(nrow(sim_res), 100)
})

test_that('Test parasite = vivax sets parameters$parasite to vivax', {
  vivax_parameters <- get_parameters(parasite = "vivax")
  expect_identical(vivax_parameters$parasite, "vivax")
})

test_that('Test difference between falciparum and vivax parameter lists', {
  falciparum_parameters <- get_parameters(parasite = "falciparum")
  vivax_parameters <- get_parameters(parasite = "vivax")
  
  expect_true(all(names(falciparum_parameters)[!names(falciparum_parameters) %in% names(vivax_parameters)] %in%
                    c("du","rvm","rva","rb","b0","b1","ib0","kb","theta0","theta1","kv","fv0","av","gammav","iv0","fd0","ad","gammad","d1","id0","kd","ub","uv","gamma1","pvm","init_ivm","init_ib","init_iva")))
  expect_true(all(names(vivax_parameters[!names(vivax_parameters) %in% names(falciparum_parameters)]) %in%
                    c("du_max","du_min","ku","au50","b","phi0lm","phi1lm","ic0lm","kclm","d0","ca","init_idm","f","gammal")))
})

## Test subpatent progression functions
test_that('Test subpatent duration function works ', {
  parameters <- get_parameters(overrides = list(s_proportion = 0,
                                                d_proportion = 0,
                                                a_proportion = 0,
                                                u_proportion = 1,
                                                t_proportion = 0))
  vivax_parameters <- get_parameters(parasite = "vivax",
                                     overrides = list(s_proportion = 0,
                                                      d_proportion = 0,
                                                      a_proportion = 0,
                                                      u_proportion = 1,
                                                      t_proportion = 0))
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
  variables$id <- individual::DoubleVariable$new(rep(1E7,100))
  variables$idm <- individual::DoubleVariable$new(rep(0,100))
  expect_equal(object = subpatent_duration(parameters = vivax_parameters, variables = variables, index = index),
               expected = rep(vivax_parameters$du_min,100),
               tolerance = 1E-2)
})

test_that('Test idm rendered creates output', {
  vivax_parameters <- get_parameters(parasite = "vivax")
  sim_res <- run_simulation(timesteps = 200, parameters = vivax_parameters)
  expect_vector(object = sim_res$idm_mean, size = 200) 
  
  expect_vector(sim_res$n_detect_730_3650)
  expect_vector(sim_res$p_detect_730_3650)
  expect_false(any(is.na(sim_res$n_detect_730_3650)))
  expect_false(any(is.na(sim_res$p_detect_730_3650)))
  
  
})


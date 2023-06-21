test_that('Test falciparum switch produces same', {
  parameters <- get_parameters()
  parameters_fal <- get_parameters(parasite = "falciparum")
  expect_equal(parameters, parameters_fal)
})

test_that('Test vivax model runs', {
  vivax_parameters <- get_parameters(parasite = "vivax_tmp")
  vivax_parameters$parasite <- "vivax"
  sim_res <- run_simulation(timesteps = 100, parameters = vivax_parameters)
  expect_equal(nrow(sim_res), 100)
})

test_that('Test parasite = vivax sets parameters$parasite to vivax', {
  vivax_parameters <- get_parameters(parasite = "vivax")
  expect_equal(vivax_parameters$parasite, "vivax")
})

test_that('Test difference between falciparum and vivax parameter lists', {
  falciparum_parameters <- get_parameters(parasite = "falciparum")
  vivax_parameters <- get_parameters(parasite = "vivax")
  
  expect_true(all(names(falciparum_parameters)[!names(falciparum_parameters) %in% names(vivax_parameters)] %in%
               c("du","rvm","rva","rb","b0","b1","ib0","kb","theta0","theta1","kv","fv0","av","gammav","iv0","fd0","ad","gammad","ub","uv","gamma1","pvm","init_ivm","init_ib","init_iva")))
  expect_true(all(names(vivax_parameters[!names(vivax_parameters) %in% names(falciparum_parameters)]) %in%
                    c("du_max","du_min","ku","au50","b","d0","ca","init_idm","f","gammal")))
})

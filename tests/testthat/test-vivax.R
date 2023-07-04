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

test_that('Test age structure should not change vivax infectivity', {
  
  # Set all individuals to asymptomatic
  # And ID immunity to not 0 (ID impacts age-specific asymptomatic infectivity)
  parameters <- get_parameters(overrides = list(s_proportion = 0,
                                                d_proportion = 0,
                                                a_proportion = 1,
                                                u_proportion = 0,
                                                t_proportion = 0,
                                                init_id  = 0.5))
  
  vivax_parameters <- get_parameters(parasite = "vivax",
                                     overrides = list(s_proportion = 0,
                                                      d_proportion = 0,
                                                      a_proportion = 1,
                                                      u_proportion = 0,
                                                      t_proportion = 0,
                                                      init_id  = 0.5))
  
  # Generate different age structure
  year <- 365
  ages <- round(c(0.083333, 1, 5, 10, 15, 20, 25, 30, 35, 40, 45,
                  50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 200) * year)
  
  deathrates <- rep(.1, length(ages)) / 365
  
  dem_parameters <- set_demography(
    parameters,
    agegroups = ages,
    timesteps = 0,
    deathrates = matrix(deathrates, nrow = 1)
  )
  
  vivax_dem_parameters <- set_demography(
    vivax_parameters,
    agegroups = ages,
    timesteps = 0,
    deathrates = matrix(deathrates, nrow = 1)
  )
  
  # vivax asymptomatic infectivity should equal ca
  expect_true(all(
    c(create_variables(vivax_parameters)$infectivity$get_values(),
      create_variables(vivax_dem_parameters)$infectivity$get_values())==vivax_parameters$ca))

  # falciparum asymptomatic infectivity should not equal ca
  expect_false(any(
    c(create_variables(parameters)$infectivity$get_values(),
      create_variables(dem_parameters)$infectivity$get_values())==vivax_parameters$ca))

  # falciparum asymptomatic infectivity should change with age structure
  expect_false(any(
    c(create_variables(parameters)$infectivity$get_values() == create_variables(dem_parameters)$infectivity$get_values())))
})

test_that('vivax schedule_infections correctly schedules new infections', {
  parameters <- get_parameters(list(human_population = 20), parasite = "vivax")
  variables <- create_variables(parameters)

  infections <- individual::Bitset$new(20)$insert(1:20)
  clinical_infections <- individual::Bitset$new(20)$insert(5:15)
  treated <- individual::Bitset$new(20)$insert(7:12)
  
  infection_mock <- mockery::mock()
  
  mockery::stub(schedule_infections, 'update_infection', infection_mock)
  
  schedule_infections(
    variables,
    clinical_infections,
    treated,
    infections,
    parameters,
    42 
  )
  
  actual_infected <- mockery::mock_args(infection_mock)[[1]][[5]]$to_vector()
  actual_asymp_infected <- mockery::mock_args(infection_mock)[[2]][[5]]$to_vector()
  
  expect_equal(
    actual_infected,
    c(5, 6, 13, 14, 15)
  )

  expect_equal(
    actual_asymp_infected,
    c(1, 2, 3, 4, 16, 17, 18, 19, 20)
  )
})

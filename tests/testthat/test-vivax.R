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
  parameters <- get_parameters(overrides =   list(s_proportion = 0,
                                                  d_proportion = 0,
                                                  a_proportion = 1,
                                                  u_proportion = 0,
                                                  t_proportion = 0,
                                                  init_id  = 0.5)
  )
  
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
  
  # Turn each of these to use vivax model: 
  vivax_parameters <- set_parasite(parameters = parameters, parasite = "vivax")
  vivax_dem_parameters <- set_parasite(parameters = dem_parameters, parasite = "vivax")

  # Run simulations for a single timestep
  sim <- run_simulation(timesteps = 1, parameters = parameters)
  sim_d <- run_simulation(timesteps = 1, parameters = dem_parameters)
  sim_v <- run_simulation(timesteps = 1, parameters = vivax_parameters)
  sim_v_d <- run_simulation(timesteps = 1, parameters = vivax_dem_parameters)
  
  # Compare infectivity
  if(sim$infectivity == sim_d$infectivity){
    stop("Something has gone wrong with the falciparum code, likely the part that 
    calculates asymptomatic infectivity")}
  
  expect_equal(sim_v$infectivity, sim_v_d$infectivity)
})

test_that('Test falciparum model runs', {
  parameters <- get_parameters()
  sim_res <- run_simulation(timesteps = 100, parameters = parameters)
  expect_equal(nrow(sim_res), 100)
})

test_that('Test model vivax runs', {
  parameters <- get_parameters()
  vivax_parameters <- set_parasite(parameters = parameters, parasite = "vivax")
  sim_res <- run_simulation(timesteps = 100, parameters = parameters)
  expect_equal(nrow(sim_res), 100)
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

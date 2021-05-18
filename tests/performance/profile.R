year <- 365
sim_length <- 1 * year
human_population <- 1e5
eir <- 1e3

simparams <- malariasimulation::get_parameters(
  list(
    human_population = human_population,
    individual_mosquitoes = FALSE
  )
)
simparams <- malariasimulation::set_equilibrium(simparams, eir)

profvis::profvis({output <- malariasimulation::run_simulation(sim_length, simparams)})
# output <- malariasimulation::run_simulation(sim_length, simparams)
print('done')
library(malariasimulation)
library(malariaEquilibrium)
library(profvis)

year <- 365
sim_length <- 1 * year
human_population <- 1e3
eir <- 1e3
ft <- .5

jamie_params <- load_parameter_set()

simparams <- get_parameters(c(
  translate_jamie(remove_unused_jamie(jamie_params)),
  list(
    human_population = human_population,
    species = 'All',
    species_proportions = 1
  )
))

eq <- human_equilibrium(EIR = eir, ft = ft, p = jamie_params, age = 0:99)
simparams <- parameterise_human_equilibrium(simparams, eq)
simparams <- parameterise_mosquito_equilibrium(simparams, EIR=eir)
simparams <- set_drugs(simparams, list(AL_params, DHC_PQP_params))
simparams <- set_clinical_treatment(simparams, ft, c(1, 2), c(.5, .5))

profvis({output <- run_simulation(sim_length, simparams)})
#output <- run_simulation(sim_length, simparams)
print('done')

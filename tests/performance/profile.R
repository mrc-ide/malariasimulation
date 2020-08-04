library(malariasimulation)
library(malariaEquilibrium)
library(profvis)

year <- 365
sim_length <- 5 * year
human_population <- 60000

jamie_params <- load_parameter_set()

simparams <- get_parameters(c(
  translate_jamie(remove_unused_jamie(jamie_params)),
  list(
    human_population = human_population,
    variety_proportions = 1,
    total_M = human_population * 5,
    mosquito_limit = human_population * 10
  )
))

simparams <- set_drugs(simparams, list(AL_params, DHC_PQP_params))
simparams <- set_clinical_treatment(simparams, .5, c(1, 2), c(.5, .5))

profvis({output <- run_simulation(sim_length, simparams)})

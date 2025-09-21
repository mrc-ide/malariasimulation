#compare runs against malariasim 

devtools::load_all()
year <- 365
sim_length <- 365*10
human_population <- 1000
starting_EIR <- 150
simparams <- get_parameters(overrides = list(
  human_population = human_population, 
  prevalence_rendering_min_ages = 0, 
  prevalence_rendering_max_ages = 5 * 365, 
  individual_mosquitoes = FALSE
))


mosq_params <- gamb_params 
mosq_params$Q0 <- 0.9
mosq_params$phi_bednets <- 0.89

simparams <- set_species(simparams, species = list(mosq_params), 
                         proportions = c(1))

simparams <- set_equilibrium(parameters = simparams, 
                             init_EIR = starting_EIR)


#then add interventions 

#add bednets
bednetstimesteps <- c(1, 4, 7, 10)*365
bednetparams <- set_bednets(
  simparams, 
  timesteps = bednetstimesteps, 
  coverages = c(0.75, 0.75, 0.75, 0.75), 
  retention = 1000*365, 
  dn0 = matrix(c(0.41, 0.41, 0.41, 0.41), nrow = 4, ncol = 1), 
  rn = matrix(c(0.56, 0.56, 0.56, 0.56), nrow = 4, ncol = 1), 
  rnm = matrix(c(0.24, 0.24, 0.24, 0.24), nrow = 4, ncol = 1), 
  gamman = rep((2.64/log(2))*365, 4)
)


correlations <- get_correlation_parameters(bednetparams)
correlations$inter_round_rho('bednets', 1)

output_reps <- run_simulation_with_repetitions(timesteps = sim_length, 
                                               overrides = list(bednetparams, correlations),
                                               repetitions = 30)
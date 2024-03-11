devtools::load_all()

require(tidyverse)
year <- 365
sim_length <- 6 * year
human_population <- 1000
starting_EIR <- 100

simparams <- get_parameters(
  list(human_population = human_population, 
       endec = TRUE)
  #sim_length = sim_length, )
)

simparams <- set_equilibrium(parameters = simparams, init_EIR = starting_EIR)

steps <- (3*365):(4*365)

mosq_params <- gamb_params
#mosq_params$Q0 <- gamb_params$Q0*1/3 #when add this in, the mu goes really low.....
simparams <- set_species(simparams, species = list(mosq_params), 
                         proportions = c(1))

endec_params <- set_endectocide(parameters = simparams, timesteps = steps, coverages = rep(1, length(steps)))

output <- run_simulation(time_length(sim_length), parameters = endec_params)

ggplot(output, aes(x = timestep, y = total_M_gamb))+
  geom_line()+
  ylim(0, 10e4)
#get a drop in mv when we have Q0*cov, but when we do Q0*1/3*cov, get an increase in mv.Even if we artificially reduce Q0 by a factor of 1/3, doesn't work
#think we don't need to account for the 1 bite every 3 days because it is captured through tau1 and tau2
ggplot(output, aes(x = timestep, y = n_detect_730_3650/n_730_3650))+
  geom_line()+ 
  ylim(0.2, 1)

ggplot(output, aes(x = timestep, y = mu_gamb))+
  geom_line()

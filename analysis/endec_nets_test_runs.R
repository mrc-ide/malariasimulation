#first comparing the effects of adding the waning effect to mortality in malariasim vs the other models

devtools::load_all()

require(tidyverse)
year <- 365
sim_length <- 1 * year
human_population <- 1e5
starting_EIR <- 60



simparams <- get_parameters(overrides = list(
  human_population = human_population,
  endec = TRUE,
  prevalence_rendering_min_ages = 0,    
  
  prevalence_rendering_max_ages = 5 * 365, 
  individual_mosquitoes = FALSE))

simparams <- set_equilibrium(parameters = simparams, init_EIR = starting_EIR) ##is this the right order? Should we set to eqm later on?
mosq_params <- gamb_params
mosq_params$Q0 <- 0.9
simparams <- set_species(simparams, species = list(mosq_params), 
                         proportions = c(1))

IVM_begin1 <- 100 # d100
mda_int <- 30
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)
eff_len <- 23
steps <- (IVM_start[1]):(IVM_start[3]+eff_len)
#steps <- (IVM_start[1]): (IVM_start[2])
#coverage_in <- c(rep(0.9, eff_len), rep(0, mda_int-eff_len), 
#                 rep(0.9, eff_len), rep(0, mda_int-eff_len), 
#                 rep(0.9, eff_len+1)) #90% coverage for first 23 days of each month


#endec_on <- c(rep(IVM_start[1], eff_len), rep(0, mda_int - eff_len),
#            rep(IVM_start[2], eff_len), rep(0, mda_int - eff_len), 
#            rep(IVM_start[3], eff_len+1))

endec_on <- c(100:123, 130:153, 160:183)
endec_ts <- c(100, 130, 160)

endec_params <- set_endectocide(parameters = simparams, timesteps = steps,
                                 endec_on = endec_on, endec_ts = endec_ts) #Error in if (expected_bites > 0) { : 


#bednetstimesteps <- c(1, 4) * year # The bed nets will be distributed at the end of the first and the 4th year. 
#
#bednetparams <- set_bednets(
#  endec_params,
#  timesteps = bednetstimesteps,
#  coverages = c(.5, .5),  # Each round is distributed to 50% of the population.
#  retention = 5 * year, # Nets are kept on average 5 years
#  dn0 = matrix(c(.533, .533), nrow = 2, ncol = 1), # Matrix of death probabilities for each mosquito species over time
#  rn = matrix(c(.56, .56), nrow = 2, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
#  rnm = matrix(c(.24, .24), nrow = 2, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
#  gamman = rep(2.64 * 365, 2) # Vector of bed net half-lives for each distribution timestep
#)

output_endec <- run_simulation_with_repetitions(timesteps = sim_length, overrides = endec_params, repetitions = 30)

ggplot(output_endec, aes(x = timestep, y = mu_gamb))+
  geom_line()+
 # geom_hline(yintercept = 0.132, col = "red")+
  geom_vline(xintercept = IVM_start[1], col = "red", linetype = "dashed")+
  geom_vline(xintercept = IVM_start[2], col = "red", linetype = "dashed")+
  geom_vline(xintercept = IVM_start[3], col = "red", linetype = "dashed")+
  theme_bw()

ggplot(output_endec, aes(x = timestep, y = total_M_gamb))+
  geom_line()+
  geom_vline(xintercept = IVM_start[1], col = "red")+
  geom_vline(xintercept = IVM_start[2], col = "red")+
  geom_vline(xintercept = IVM_start[3], col = "red")

ggplot(output_endec, aes(x = timestep, y = EIR_gamb))+
  geom_line()+
  geom_hline(yintercept = 0.132, col = "red")

ggplot(output_endec, aes(x = timestep, y = n_detect_0_1825/n_0_1825))+
  geom_line()+
  ylim(0, 0.8) #malariasim slightly overestimates. 

#read in the odin models
odin_models <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/exploring_interactions/malariasim-odin/malariasim-odin-compare2.rds")
odin_models <- odin_models %>%
  filter(ivermectin_coverage == "high") %>%
  select(t, mv, EIR_tot, slide_prev0to5, Ivtot, model_type)

output_endec2 <- output_endec %>%
  mutate(slide_prev0to5 = n_detect_lm_0_1825/n_age_0_1825) %>%
  select(timestep, total_M_gamb, EIR_gamb, slide_prev0to5, Im_gamb_count) %>%
  mutate(model_type = "malariasim") %>%
  rename(t = timestep, 
         mv = total_M_gamb, 
         EIR_tot = EIR_gamb, 
         Ivtot = Im_gamb_count)

head(output_endec2)
head(odin_models)

output <- rbind(output_endec2, odin_models)

prev_plot <- ggplot(output, aes(x = t, y = slide_prev0to5, col = as.factor(model_type)))+
  geom_line(linewidth =1.5, alpha = 0.8)+
  theme_bw()+
  ylim(0, 1)


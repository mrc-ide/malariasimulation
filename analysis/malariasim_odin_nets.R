#compare runs against malariasim 

devtools::load_all()
year <- 365
sim_length <- 365*10
human_population <- 1000
starting_EIR <- 150
simparams <- get_parameters(overrides = list(
  human_population = human_population, 
  endec = TRUE, 
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

#add endectocide

IVM_begin1 <- 365*2
mda_int <- 30 
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)
eff_len <- 23
steps <- (IVM_start[1]):(IVM_start[3]+eff_len)
endec_on <- c(as.numeric(IVM_start[1]):as.numeric(IVM_start[1]+eff_len), 
              as.numeric(IVM_start[2]):as.numeric(IVM_start[2]+eff_len), 
              as.numeric(IVM_start[3]):as.numeric(IVM_start[3]+eff_len))
endec_ts <- c(IVM_start[1], IVM_start[2], IVM_start[3])

endec_params <- set_endectocide(parameters = bednetparams, timesteps = steps,
                                endec_on = endec_on, endec_ts = endec_ts)


output <- run_simulation(timesteps = sim_length, parameters = endec_params, correlations = correlations)

#output_reps <- run_simulation_with_repetitions(timesteps = sim_length, 
#                                               overrides = list(endec_params, correlations),
#                                               repetitions = 30)
#
#
#----- 2) Version with rewritten run_simulations_with_repetitions() function -----------------------

##' The updated version of run_simulation_with_repetitions() replaces the overrides argument with
##' parameters and correlations arguments, which are now fed directly to the internal run_simulations()
##' function call.

# Store the updated version of run_simulation_with_repetitions():
#Tom Brewer
updated_run_simulation_with_repetitions <- function(
    timesteps,
    repetitions,
    parameters,
    correlations,
    parallel = FALSE
) {
  if (parallel) {
    fapply <- parallel::mclapply
  } else {
    fapply <- lapply
  }
  dfs <- fapply(
    seq(repetitions),
    function(repetition) {
      df <- run_simulation(
        timesteps = timesteps,
        parameters = parameters,
        correlations = correlations
      )
      df$repetition <- repetition
      df
    }
  )
  do.call("rbind", dfs)
}

# Run the simulations using the updated_run_simulation_with_repetitions() function:
output_reps <- updated_run_simulation_with_repetitions(
  timesteps = sim_length,
  repetitions = 30,
  parameters = endec_params,
  correlations = correlations,
  parallel = FALSE
)

ggplot(output, aes(x = timestep, y = n_use_net/human_population))+
  geom_line()+
  ylim(0,1)


require(tidyverse)
ggplot(output, aes(x = timestep, y = n_detect_lm_0_1825/n_age_0_1825))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()

ggplot(output, aes(x = timestep, y = total_M_gamb))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()

ggplot(output, aes(x = timestep, y = (EIR_gamb*365)/1000))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()



saveRDS(output_reps, file = "analysis/output_reps_net_endec.rds")

output_reps <- readRDS("analysis/output_reps_net_endec.rds")

output_reps_summary <- output_reps %>% #otherwise output_reps
  group_by(timestep) %>%
  summarise(mu = mean(mu_gamb), 
            mv = mean((total_M_gamb)/1000), 
            EIR_tot = (mean(EIR_gamb)*365)/1000, 
            slide_prev0to5 = mean(n_detect_lm_0_1825/n_age_0_1825))%>%
  
  mutate(model = "malariasim") %>%
  rename(t = timestep) %>%
  select(t, mv, EIR_tot, slide_prev0to5, model, mu)

#output_summary <- output %>% #otherwise output_reps
#  group_by(timestep) %>%
#  summarise(mu = mean(mu_gamb), 
#            mv = mean(total_M_gamb), 
#            EIR_tot = (mean(EIR_gamb)*365)/1000, 
#            slide_prev0to5 = mean(n_detect_lm_0_1825/n_age_0_1825))%>%
#  
#  mutate(model = "malariasim") %>%
#  rename(t = timestep) %>%
#  select(t, mv, EIR_tot, slide_prev0to5, model)

#output %>%
#  group_by(timestep) %>%
#  summarise(net_use = mean(n_use_net/human_population)) %>%
#  ggplot()+
#  aes(x = timestep, y = net_use)+
#  geom_line()

ggplot(output_reps_summary, aes(x = t, y = slide_prev0to5))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()

ggplot(output_reps_summary, aes(x = t, y = mv))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()

ggplot(output_reps_summary, aes(x = t, y = EIR_tot))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()


#made in ivRmectin...ivm_nets_malariasim.R
odin_models <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/exploring_interactions/malariasim-odin/hazards_malariasim_odin_nets.rds")

odin_models <- odin_models %>%
  select(t, mv, EIR_tot, slide_prev0to5, model, mu)

all_models <- rbind(odin_models, output_reps_summary)


all_models %>%
  filter(model != "malariasim-odin-exp-decay") %>%
  ggplot()+
  aes(x = t, y = EIR_tot, col = as.factor(model))+
  geom_line(size = 1)+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "black", linetype = "dashed")+
  theme_minimal()

ggplot(all_models, aes(x = t, y = EIR_tot, col = as.factor(model)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()

ggplot(all_models, aes(x = t, y = slide_prev0to5, col = as.factor(model)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()+
  ylim(0,1)

ggplot(all_models, aes(x = t, y = mv, col = as.factor(model)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()

ggplot(all_models, aes(x = t, y = mu, col = as.factor(model)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()+
  ylim(0,1)

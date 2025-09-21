#compare runs against malariasim 
require(tidyverse)


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
  repetitions = 5,
  parameters = bednetparams,
  correlations = correlations,
  parallel = FALSE
)

write_rds(output_reps, file = "output_reps_nets_only.rds")

output_reps <- readRDS("output_reps_nets_only.rds")


ggplot(output_reps, aes(x = timestep, y = mu_gamb, col = as.factor(repetition)))+
  geom_line()


output_reps_summary <- output_reps %>% #otherwise output_reps
  group_by(timestep) %>%
  summarise(mu = mean(mu_gamb), 
            mv = mean(total_M_gamb),
            #mv = mean(total_M_gamb)*omega/1000, 
           # mv = mean(total_M_gamb/total_M_gamb[1]), 
            EIR_tot = (mean(EIR_gamb)*365)/1000, 
            prev = mean(n_detect_lm_0_1825/n_age_0_1825), 
            FOIm = mean(FOIM_gamb))
  
  mutate(model = "malariasim") %>%
  rename(t = timestep) %>%
  select(t, mv, EIR_tot, prev, model, mu, FOIm)


ggplot(output_reps_summary, aes(x = t, y = FOIm))+
  geom_line()+
  ylim(0,0.025)

ggplot(output_reps_summary, aes(x = t, y = mv/mv[1]))+
  geom_line()

ggplot(output_reps_summary, aes(x = t, y = mv))+
  geom_line()

#don't take mean#
output_reps2 <- output_reps %>% #otherwise output_reps
  group_by(timestep) %>%
  summarise(mu = mu_gamb, 
            mv = total_M_gamb,
            #mv = mean(total_M_gamb)*omega/1000, 
            # mv = mean(total_M_gamb/total_M_gamb[1]), 
            EIR_tot = (EIR_gamb*365)/1000, 
            prev = n_detect_lm_0_1825/n_age_0_1825, 
            FOIm = FOIM_gamb, 
            repetition = repetition) %>%
  mutate(model = "malariasim") %>%
  rename(t = timestep) %>%
  select(t, mv, EIR_tot, prev, model, mu, FOIm, repetition) %>%
  mutate(repetition = case_when(repetition == 1 ~ "malariasim_rep1", 
                                repetition == 2 ~ "malariasim_rep2", 
                                repetition == 3 ~ "malariasim_rep3", 
                                repetition == 4 ~ "malariasim_rep4", 
                                repetition == 5 ~ "malariasim_rep5"))


#made in ICDMM nets_runs.R
odin_nets <- readRDS("C:/Users/nc1115/Documents/github/deterministic-malaria-model/net_out_odin.rds")
omega_odin <- unique(odin_nets$omega)
odin_nets <- odin_nets %>%
  mutate(model = "odin") %>%
  select(t, mv, EIR_tot, prev, model, mu, FOIv) %>%
  rename(FOIm = FOIv) %>%
  mutate(mv = mv*1000/omega_odin, 
         repetition = "Odin")



all_models <- rbind(odin_nets, output_reps2)

foi_plot <- ggplot(all_models, aes(x = t, y = FOIm, col = as.factor(repetition)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "black", linetype = "dashed")+
  theme_minimal()+
  ylim(0, 0.01)




mu_plot <- ggplot(all_models, aes(x = t, y = mu, col = as.factor(repetition)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "black", linetype = "dashed")+
  theme_minimal()

mv_plot <- ggplot(all_models, aes(x = t, y = mv, col = as.factor(repetition)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "black", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "black", linetype = "dashed")+
  theme_minimal()


eir_plot <- ggplot(all_models, aes(x = t, y = EIR_tot, col = as.factor(repetition)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
 #geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()


prev_plot <- ggplot(all_models, aes(x = t, y = prev, col = as.factor(repetition)))+
  geom_line()+
  geom_vline(xintercept = 365*1, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*4, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  geom_vline(xintercept = 365*7, col = "blue", linetype = "dashed")+
  #geom_vline(xintercept = 365*2, col = "red", linetype = "dashed")+
  theme_minimal()



cowplot::plot_grid(prev_plot, mv_plot, eir_plot, foi_plot, mu_plot)

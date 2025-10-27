devtools::load_all()
require(tidyverse)
#malariasim with baseline scenario, ITN and endectocides

#baseline scenario####

year <- 365
sim_length <- 365*15
human_population <- 1e5
starting_EIR <- 115


# Set the age ranges (in days)
age_min <- seq(0, 80, 5) * 365
age_max <- seq(5, 85, 5) * 365


simparams <- get_parameters(overrides = list(
  human_population = human_population, 
  endec = TRUE, 
  prevalence_rendering_min_ages = 0, 
  prevalence_rendering_max_ages = 5 * 365, 
  clinical_incidence_rendering_min_ages = 0, 
  clinical_incidence_rendering_max_ages = 5*365,
  individual_mosquitoes = FALSE, 
  age_group_rendering_min_ages = age_min, 
  
  age_group_rendering_max_ages = age_max
))

mosq_params <- gamb_params 
mosq_params$Q0 <- 0.95
mosq_params$phi_bednets <- 0.95 #phi-bed from Bijagos

simparams <- set_species(simparams, species = list(mosq_params), 
                         proportions = c(1))

simparams <- set_equilibrium(parameters = simparams, 
                             init_EIR = starting_EIR)

baseline_sim <- run_simulation(timesteps = sim_length,simparams)
#saveRDS(baseline_sim, file = "analysis/chapter-int-ecology/run_baseline_sim.rds")

odin_baseline <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/exploring_interactions/MIM_poster/chapter_plots/df_baseline.rds")
odin_high_eir <- odin_baseline %>%
  filter(ref == 2)


prev_compare <- ggplot(baseline_sim, aes(x = timestep, y = (n_detect_lm_0_1825/n_age_0_1825)*100))+
  geom_line()+
  ylim(0,75)+
  geom_line(data = odin_high_eir, aes(x = t, y = slide_prev0to5*100), col = "red")

inc_compare <- ggplot(baseline_sim, aes(x = timestep, y = (n_inc_clinical_0_1825/n_age_0_1825)*1000))+
  geom_line()+
  geom_line(data = odin_high_eir, aes(x = t, y = clin_inc0to5*1000), col = "red")+
  ylab("Clinical incidence per 1000 children aged 0-5")+
  ylim(0,10)

cowplot::plot_grid(prev_compare, inc_compare)


#gplot(baseline_sim, aes(x = timestep, y = (n_inc_clinical_0_1825/n_age_0_1825)*1e5))+
# geom_line()+
# #ylim(0,100)+
# geom_line(data = odin_high_eir, aes(x = t, y = clin_inc0to5*1000), col = "red")


##vector control####
#ITNs: setup
itn_on <- 365*5

net_seq <- seq(itn_on, sim_length, by = 3*365)
#bednet_years <- seq(3, 10, by = 3)
bednet_timesteps <- net_seq
itn_distr <- length(net_seq)

retention_net <- c(5, 2, 0.5)*365 #long, medium and short

#ITN efficacy data
net_eff_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_only_nets.csv")

res_55 <- net_eff_df %>%
  filter(resistance == 0.55)


#run with long retention
bednet_params_long <- set_bednets(simparams, 
                                  timesteps = bednet_timesteps, 
                                  retention =retention_net[1], 
                                  coverages = rep(0.5, itn_distr), #realistic for Africa
                                  #get correct dn0, rn0 etc from Ellie's work
                                  dn0 = matrix(rep(res_55$dn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of death probabilities for each mosquito species over time
                                  rn = matrix(rep(res_55$rn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                                  rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                                  gamman = rep((res_55$gamman_med/log(2)) * 365, itn_distr)) # Vector of bed net half-lives for each distribution timestep)

correlations_long <- get_correlation_parameters(bednet_params_long)
correlations_long$inter_round_rho('bednets', 1)

bednet_params_medium <- set_bednets(simparams, 
                                    timesteps = bednet_timesteps, 
                                    retention =retention_net[2], 
                                    coverages = rep(0.5, itn_distr), #realistic for Africa
                                    #get correct dn0, rn0 etc from Ellie's work
                                    dn0 = matrix(rep(res_55$dn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of death probabilities for each mosquito species over time
                                    rn = matrix(rep(res_55$rn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                                    rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                                    gamman = rep((res_55$gamman_med/log(2)) * 365, itn_distr)) # Vector of bed net half-lives for each distribution timestep)

correlations_medium <- get_correlation_parameters(bednet_params_medium)
correlations_medium$inter_round_rho('bednets', 1)

bednet_params_short <- set_bednets(simparams, 
                                    timesteps = bednet_timesteps, 
                                    retention =retention_net[3], 
                                    coverages = rep(0.5, itn_distr), #realistic for Africa
                                    #get correct dn0, rn0 etc from Ellie's work
                                    dn0 = matrix(rep(res_55$dn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of death probabilities for each mosquito species over time
                                    rn = matrix(rep(res_55$rn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                                    rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                                    gamman = rep((res_55$gamman_med/log(2)) * 365, itn_distr)) # Vector of bed net half-lives for each distribution timestep)

correlations_short <- get_correlation_parameters(bednet_params_short)
correlations_short$inter_round_rho('bednets', 1)

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



run_bednets_long <- run_simulation(parameters = bednet_params_long, 
                                                            timesteps = sim_length, 
                                                            #repetitions = 5,
                                                            correlations = correlations_long)

saveRDS(run_bednets_long, file = "analysis/chapter-int-ecology/run_bednets_long.rds")

run_bednets_medium <- run_simulation(parameters = bednet_params_medium, 
                                   timesteps = sim_length, 
                                   #repetitions = 5,
                                   correlations = correlations_medium)

saveRDS(run_bednets_medium, file = "analysis/chapter-int-ecology/run_bednets_medium.rds")

run_bednets_short <- run_simulation(parameters = bednet_params_short, 
                                     timesteps = sim_length, 
                                     #repetitions = 5,
                                     correlations = correlations_short)

saveRDS(run_bednets_short, file = "analysis/chapter-int-ecology/run_bednets_short.rds")




ggplot(run_bednets_long, aes(x = timestep/365, y = (n_detect_pcr_0_1825/n_age_0_1825)*100))+
  geom_line()+
  ylim(0, 100)+
  geom_vline(xintercept = 3, col = "red", lty = "dashed")+
  geom_vline(xintercept = 6, col = "red", lty = "dashed")+
  geom_vline(xintercept = 9, col = "red", lty = "dashed")

#short retention time
bednet_params_short <-  set_bednets(simparams, 
                                    timesteps = bednet_timesteps, 
                                    retention =retention_net[3], 
                                    coverages = rep(0.5, itn_distr),
                                    dn0 = matrix(rep(res_55$dn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of death probabilities for each mosquito species over time
                                    rn = matrix(rep(res_55$rn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                                    rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                                    gamman = rep((res_55$gamman_med/log(2)) * 365, itn_distr)) # Vector of bed net half-lives for each distribution timestep)

correlations_short <- get_correlation_parameters(bednet_params_short)
correlations_short$inter_round_rho('bednets', 1)

#medium retention 
bednet_params_medium <- set_bednets(simparams, 
                                    timesteps = bednet_timesteps, 
                                    retention =retention_net[2], 
                                    coverages = rep(0.5, itn_distr),
                                    dn0 = matrix(rep(res_55$dn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of death probabilities for each mosquito species over time
                                    rn = matrix(rep(res_55$rn0_med, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                                    rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                                    gamman = rep((res_55$gamman_med/log(2)) * 365, itn_distr)) # Vector of bed net half-lives for each distribution timestep)

correlations_medium <- get_correlation_parameters(bednet_params_medium)
correlations_medium$inter_round_rho('bednets', 1)
#endectocide####

#endectocide MDA happens early after bednet campaign
mda_int <- 30
eff_len <- 23
early_IVM <- 180
early_IVM_begin1 <- bednet_timesteps[2]+early_IVM # 6 months into bednet campaign
early_IVM_start <- c(early_IVM_begin1, early_IVM_begin1+mda_int, early_IVM_begin1+mda_int+mda_int)
early_eff_len <- 23
early_steps <- (early_IVM_start[1]):(early_IVM_start[3]+eff_len)
early_endec_on <- c(as.numeric(early_IVM_start[1]):as.numeric(early_IVM_start[1]+eff_len), 
              as.numeric(early_IVM_start[2]):as.numeric(early_IVM_start[2]+eff_len), 
              as.numeric(early_IVM_start[3]):as.numeric(early_IVM_start[3]+eff_len))
early_endec_ts <- c(early_IVM_start[1], early_IVM_start[2], early_IVM_start[3])

#endectocide MDA happens late after bednet campaign
mda_int <- 30
late_IVM <- 2*365 #2.y after
late_IVM_begin1 <- bednet_timesteps[2]+late_IVM #
late_IVM_start <- c(late_IVM_begin1, late_IVM_begin1+mda_int, late_IVM_begin1+mda_int+mda_int)
late_eff_len <- 23
late_steps <- (late_IVM_start[1]):(late_IVM_start[3]+eff_len)
late_endec_on <- c(as.numeric(late_IVM_start[1]):as.numeric(late_IVM_start[1]+eff_len), 
                    as.numeric(late_IVM_start[2]):as.numeric(late_IVM_start[2]+eff_len), 
                    as.numeric(late_IVM_start[3]):as.numeric(late_IVM_start[3]+eff_len))
late_endec_ts <- c(late_IVM_start[1], late_IVM_start[2], late_IVM_start[3])


#early MDA with long net retention
endec_params_early_endec_long_bednet <- set_endectocide(bednet_params_long, 
                                                        timesteps = early_steps,
                                                        endec_on = early_endec_on, 
                                                        endec_ts = early_endec_ts)

#early MDA with medium net retention 
endec_params_early_endec_medium_bednet <- set_endectocide(bednet_params_medium, 
                                                   timesteps = early_steps,
                                                   endec_on = early_endec_on, 
                                                   endec_ts = early_endec_ts)

#early MDA with short net retention 
endec_params_early_endec_short_bednet <- set_endectocide(bednet_params_short, 
                                                  timesteps = early_steps,
                                                  endec_on = early_endec_on, 
                                                  endec_ts = early_endec_ts)


#late MDA with long net retention 
endec_params_late_endec_long_bednet <- set_endectocide(bednet_params_long, 
                                                       timesteps = late_steps,
                                                       endec_on = late_endec_on, 
                                                       endec_ts = late_endec_ts)

#late MDA with medium net retention
endec_params_late_endec_medium_bednet <- set_endectocide(bednet_params_medium, 
                                                         timesteps = late_steps,
                                                         endec_on = late_endec_on, 
                                                         endec_ts = late_endec_ts)

#late MDA with short net retention
endec_params_late_endec_short_bednet <- set_endectocide(bednet_params_short, 
                                                        timesteps = late_steps,
                                                        endec_on = late_endec_on, 
                                                        endec_ts = late_endec_ts)

#run model####
run_endec_early_long_bednet <- run_simulation(parameters = endec_params_early_endec_long_bednet, 
                                                                       timesteps = sim_length, 
                                                                       #repetitions = 5,
                                                                       correlations = correlations_long#, 
                                                                       #parallel = FALSE
                                              )

saveRDS(run_endec_early_long_bednet, file = "analysis/chapter-int-ecology/run_endec_early_long_bednet.rds")

run_endec_early_medium_bednet <- run_simulation(parameters = endec_params_early_endec_medium_bednet, 
                                                                         timesteps = sim_length, 
                                                                         #repetitions = 5,
                                                                         correlations = correlations_medium #, 
                                                                         #parallel = FALSE
                                                                         ) 

saveRDS(run_endec_early_medium_bednet, file = "analysis/chapter-int-ecology/run_endec_early_medium_bednet.rds")

run_endec_early_short_bednet <- run_simulation(parameters = endec_params_early_endec_short_bednet, 
                                                                        timesteps = sim_length, 
                                                                        #repetitions = 5,
                                                                        correlations = correlations_short #, 
                                                                        #parallel = FALSE
                                                                        )

saveRDS(run_endec_early_short_bednet, file = "analysis/chapter-int-ecology/run_endec_early_short_bednet.rds")

run_endec_late_long_bednet <- run_simulation(parameters = endec_params_late_endec_long_bednet, 
                                             timesteps = sim_length, 
                                             #repetitions = 5,
                                             correlations = correlations_long #, 
                                             #parallel = FALSE
                                             )
saveRDS(run_endec_late_long_bednet, file = "analysis/chapter-int-ecology/run_endec_late_long_bednet.rds")

run_endec_late_medium_bednet <- run_simulation(parameters = endec_params_late_endec_medium_bednet, 
                                               timesteps = sim_length, 
                                               #repetitions = 5,
                                               correlations = correlations_medium #, 
                                               #parallel = FALSE
                                               )
saveRDS(run_endec_late_medium_bednet, file = "analysis/chapter-int-ecology/run_endec_late_medium_bednet.rds")

run_endec_late_short_bednet <- run_simulation(parameters=endec_params_late_endec_short_bednet, 
                                              timesteps = sim_length, 
                                              #repetitions = 5,
                                              correlations = correlations_short #, 
                                              #parallel = FALSE
                                              )

saveRDS(run_endec_late_short_bednet, file = "analysis/chapter-int-ecology/run_endec_late_short_bednet.rds")

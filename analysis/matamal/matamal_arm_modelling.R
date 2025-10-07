devtools::load_all()
require(tidyverse)
#malariasim with baseline scenario, ITN and endectocides

#baseline scenario####

#matamal trial modelling

year <- 365
sim_length <- 365*10
human_population <- 1e5
starting_EIR <- 150


# Set the age ranges (in days)
age_min <- seq(0, 80, 5) * 365
age_max <- seq(5, 85, 5) * 365

seasonality_params <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/1.data-cleaning/output/seasonality_params.rds")

simparams <- get_parameters(overrides = list(
  human_population = human_population, 
  endec = TRUE, 
  bednets = TRUE,
  prevalence_rendering_min_ages = 0, 
  prevalence_rendering_max_ages = 85 * 365, #all age prevalence
  individual_mosquitoes = FALSE, 
  age_group_rendering_min_ages = age_min, 
  age_group_rendering_max_ages = age_max,
  model_seasonality = TRUE,
  g0 = seasonality_params$g0,
  g = c(seasonality_params$g1, seasonality_params$g2, seasonality_params$g3),
  h = c(seasonality_params$h1, seasonality_params$h2, seasonality_params$h3)
  
))

mosq_params <- gamb_params 

#setting parameters
df_setting <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_odin_inputs.rds")
df_phi <- read.csv("analysis/matamal/phi_estimates_bijagos.csv", header = TRUE) %>%
  filter(scenario == "Data")


mosq_params$Q0 <- df_setting$hbi_median[1]
mosq_params$phi_bednets <- round(df_phi$mean,3) #phi-bed from Bijagos

simparams <- set_species(simparams, species = list(mosq_params), 
                         proportions = c(1))

simparams <- set_equilibrium(parameters = simparams, 
                             init_EIR = starting_EIR)



ggplot(baseline_sim, aes(x = timestep, y = (n_detect_pcr_0_1825/n_age_0_1825)*100))+
  geom_line()+
  ylim(0,100)

##vector control####
#ITNs: setup

bednet_years <- seq(3, 10, by = 3)
bednet_timesteps <- bednet_years*year
itn_distr <- length(bednet_years)

retention_net <- c(5, 2, 0.5) #long, medium and short


#run with long retention
bednet_params_long <- set_bednets(simparams, 
                                  timesteps = bednet_timesteps, 
                                  retention =retention_net[1], 
                                  coverages = rep(0.8, itn_distr),
                                  dn0 = matrix(rep(0.533, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of death probabilities for each mosquito species over time
                                  rn = matrix(rep(0.56, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                                  rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                                  gamman = rep(2.64 * 365, itn_distr)) # Vector of bed net half-lives for each distribution timestep)

correlations_long <- get_correlation_parameters(bednet_params_long)
correlations_long$inter_round_rho('bednets', 1)



run_bednets_long <- run_simulation(bednet_params_long, timesteps = sim_length, 
                                   correlations = correlations_long)

saveRDS(run_bednets_long, file = "analysis/chapter-int-ecology/run_bednets_long.rds")



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
                                    coverages = rep(0.8, itn_distr),
                                    dn0 = matrix(rep(0.533, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of death probabilities for each mosquito species over time
                                    rn = matrix(rep(0.56, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                                    rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                                    gamman = rep(2.64 * 365, itn_distr)) # Vector of bed net half-lives for each distribution timestep)

correlations_short <- get_correlation_parameters(bednet_params_short)
correlations_short$inter_round_rho('bednets', 1)




#medium retention 
bednet_params_medium <- set_bednets(simparams, 
                                    timesteps = bednet_timesteps, 
                                    retention =retention_net[2], 
                                    coverages = rep(0.8, itn_distr),
                                    dn0 = matrix(rep(0.533, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of death probabilities for each mosquito species over time
                                    rn = matrix(rep(0.56, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                                    rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                                    gamman = rep(2.64 * 365, itn_distr)) # Vector of bed net half-lives for each distribution timestep)

correlations_medium <- get_correlation_parameters(bednet_params_medium)
correlations_medium$inter_round_rho('bednets', 1)
#endectocide####

#endectocide MDA happens early after bednet campaign
mda_int <- 30
eff_len <- 23
early_IVM <- 180
early_IVM_begin1 <- bednet_timesteps[1]+early_IVM # 6 months into bednet campaign
early_IVM_start <- c(early_IVM_begin1, early_IVM_begin1+mda_int, early_IVM_begin1+mda_int+mda_int)
early_eff_len <- 23
early_steps <- (early_IVM_start[1]):(early_IVM_start[3]+eff_len)
early_endec_on <- c(as.numeric(early_IVM_start[1]):as.numeric(early_IVM_start[1]+eff_len), 
                    as.numeric(early_IVM_start[2]):as.numeric(early_IVM_start[2]+eff_len), 
                    as.numeric(early_IVM_start[3]):as.numeric(early_IVM_start[3]+eff_len))
early_endec_ts <- c(early_IVM_start[1], early_IVM_start[2], early_IVM_start[3])

#endectocide MDA happens late after bednet campaign
mda_int <- 30
late_IVM <- 2*365 #2y after
late_IVM_begin1 <- bednet_timesteps[1]+late_IVM # 6 months into bednet campaign
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
                                                         endec_on = late_endec_on, 
                                                         endec_ts = late_endec_ts)


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
run_endec_early_long_bednet <- run_simulation(endec_params_early_endec_long_bednet, 
                                              timesteps = sim_length, 
                                              correlations = correlations_long)

saveRDS(run_endec_early_long_bednet, file = "analysis/chapter-int-ecology/run_endec_early_long_bednet.rds")

run_endec_early_medium_bednet <- run_simulation(endec_params_early_endec_medium_bednet, 
                                                timesteps = sim_length, 
                                                correlations = correlations_medium) 
saveRDS(run_endec_early_medium_bednet, file = "analysis/chapter-int-ecology/run_endec_early_medium_bednet.rds")

run_endec_early_short_bednet <- run_simulation(endec_params_early_endec_short_bednet, 
                                               timesteps = sim_length, 
                                               correlations = correlations_short)

saveRDS(run_endec_early_short_bednet, file = "analysis/chapter-int-ecology/run_endec_early_short_bednet.rds")

run_endec_late_long_bednet <- run_simulation(endec_params_late_endec_long_bednet, 
                                             timesteps = sim_length, 
                                             correlations = correlations_long)
saveRDS(run_endec_late_long_bednet, file = "analysis/chapter-int-ecology/run_endec_late_long_bednet.rds")

run_endec_late_medium_bednet <- run_simulation(endec_params_late_endec_medium_bednet, 
                                               timesteps = sim_length, 
                                               correlations = correlations_medium)
saveRDS(run_endec_late_medium_bednet, file = "analysis/chapter-int-ecology/run_endec_late_medium_bednet.rds")

run_endec_late_short_bednet <- run_simulation(endec_params_late_endec_short_bednet, 
                                              timesteps = sim_length, 
                                              correlations = correlations_short)

saveRDS(run_endec_late_short_bednet, file = "analysis/chapter-int-ecology/run_endec_late_short_bednet.rds")

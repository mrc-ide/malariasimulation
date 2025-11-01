#sensitivity analysis on dynamics for median ivermectin coverage and Q0

#parameter set - for these, endec_mu is 0.15 and wane is 0.01
df_scenario <- readRDS("analysis/matamal/output/df_test_scenario.rds")

devtools::load_all()
require(tidyverse)
require(cali)
#malariasim with baseline scenario, ITN and endectocides

#baseline scenario####

#matamal trial modelling

year <- 365
sim_length <- 365*13 #From 2012 
human_population <- 1e5
#starting_EIR <- 10


# Set the age ranges (in days)
age_min <- seq(0, 75, 5) * 365
age_max <- seq(5, 80, 5) * 365





seasonality_params <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/1.data-cleaning/output/seasonality_params.rds")



site_demo <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/1.data-cleaning/output/site_file_demography.rds")
site_demo_2012 <- site_demo %>%
  filter(year >= 2012)

dim(site_demo_2012)

add_demography <- function(p, demography){
  
  # Age group upper
  ages <- round(unique(demography$age_upper) * 365)
  timesteps <-  (unique(demography$year) - 2012)*365 #2012 or whatever year you are starting your simulation
  deathrates <- demography$adj_mort_rates / 365
  deathrates_matrix <- matrix(deathrates, nrow = length(timesteps), byrow = TRUE)
  # Add parameters
  p <- malariasimulation::set_demography(
    parameters = p,
    agegroups = ages,
    timesteps = timesteps,
    deathrates = deathrates_matrix
  )
  
  return(p)
}


mosq_params <- gamb_params
sim_list <- list()
df_setting <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_odin_inputs.rds")

timestep_from_2012 <- function(year, month, day) {
  # Days in each month (non-leap year)
  month_lengths <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  # Years since 2012 (each year = 365 days)
  years_since <- year - 2012
  days_from_years <- years_since * 365
  
  # Days from previous months in the same year
  if (month > 1) {
    days_from_months <- sum(month_lengths[1:(month - 1)])
  } else {
    days_from_months <- 0
  }
  
  # Days within the current month (start at day 1 → offset 0)
  days_from_days <- day - 1
  
  # Total timestep (timestep = 1 corresponds to 1 Jan 2012)
  timestep <- 1 + days_from_years + days_from_months + days_from_days
  return(timestep)
}

#test with one row 
#df_scenario <- df_scenario[1,]

for (i in seq_len(nrow(df_scenario))) {
  Q0_i <- df_scenario$Q0[i]
  phi_bednets_i <- df_scenario$bites_Bed[i]
  retention_i <- df_scenario$retention_net_in[i]
  endec_mu_i <- df_scenario$endec_mu[i]
  wane_endec_i <- df_scenario$wane_endec[i]
  
  # update mosquito parameters
  mosq_params_i <- mosq_params
  mosq_params_i$Q0 <- Q0_i
  mosq_params_i$phi_bednets <- phi_bednets_i
  
  #parameter list
  simparams1 <- get_parameters(overrides = list(
    human_population = human_population, 
    endec = TRUE, 
    bednets = TRUE,
    mu_endec = c(endec_mu_i, endec_mu_i, endec_mu_i), 
    wane_endec = c(wane_endec_i, wane_endec_i, wane_endec_i),
    prevalence_rendering_min_ages = c(0*365,0*365),
    prevalence_rendering_max_ages = c(5*365,80*365), #all age prevalence - and under 5s
    individual_mosquitoes = FALSE, 
    age_group_rendering_min_ages = age_min, 
    age_group_rendering_max_ages = age_max ,
    model_seasonality = TRUE,
    g0 = seasonality_params$g0,
    g = c(seasonality_params$g1, seasonality_params$g2, seasonality_params$g3),
    h = c(seasonality_params$h1, seasonality_params$h2, seasonality_params$h3)
    
  ))
  
  # set the demography 
  demo_params <- add_demography(p = simparams1, demography=site_demo_2012)

  
  # set species
  sim_params_i <- set_species(
    demo_params, species = list(mosq_params_i), 
    proportions = c(1)
  )
  
  #set drugs 
  drug_params_i <- set_drugs(sim_params_i, list(AL_params))
  
  sf_treatment <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/1.data-cleaning/output/site_file_tx_cov_AL.rds")
  names(sf_treatment) <- c("year", "tx_cov")
  sf_treatment_2012 <- sf_treatment %>%
    filter(year >= 2012)
  
  
  
  # Set treatment program for AL (drug index = 1)
  treatment_params_i <- set_clinical_treatment(
    parameters = drug_params_i,
    drug = 1, #just AL
    timesteps =  (sf_treatment_2012$year - 2012)*365 + 1, # Treatment coverage changes on day 300 and day 600
    coverages =  sf_treatment_2012$tx_cov) # The initial treatment coverage (0%) is the default. site file would multiply by prop_ACT and here we assume that to be 1  
  
  
  #set bednets
  df_nets <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_net_matamal_info.rds")
  df_nets2 <- df_nets  %>%
    filter(year >= 2012)
  
  itn_distr <- nrow(df_nets2)
  
  distr_campaign <- 6*30 # mass distributions are around June
  
  med_dn0 <- unique(df_nets2$dn0_med)
  med_rn0 <- unique(df_nets2$rn0_med)
  med_gamman <- unique(df_nets2$gamman_med)
  
  #some sort of bug here in set_bednets
  bednet_params_i <- set_bednets(treatment_params_i, 
                               timesteps = (df_nets2$year - 2012)*365 + distr_campaign, 
                               retention = retention_i, 
                               coverages = df_nets2$itn_input_distr,
                               dn0 = matrix(rep(med_dn0, itn_distr), nrow = itn_distr,ncol = 1), # Matrix of death probabilities for each mosquito species over time
                               rn = matrix(rep(med_rn0, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                               rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                               gamman = rep(med_gamman * 365, itn_distr))
  
  correlations <- get_correlation_parameters(bednet_params_i)
  correlations$inter_round_rho('bednets', 1)
  bednet_params_i$timesteps <- sim_length
  input_EIR <- 7
  simparams_i <- set_equilibrium(parameters = bednet_params_i, 
                               init_EIR = input_EIR)
  
  
  
  
  #then add the DP 
  mdaparams_i <- simparams_i
  mdaparams_i <- set_drugs(mdaparams_i, list(DHA_PQP_params))
  
  DP_july_21 <- timestep_from_2012(2021, 7, 1)
  DP_aug_21 <- timestep_from_2012(2021, 8, 1)
  DP_sept_21 <- timestep_from_2012(2021, 9, 1)
  
  DP_july_22 <- timestep_from_2012(2022, 7, 1)
  DP_aug_22 <- timestep_from_2012(2022, 8, 1)
  DP_sept_22 <- timestep_from_2012(2022, 9, 1)
  
  #DP_mda_y_2021 <- c(y_2021+july, y_2021+aug, y_2021+sept)
  DP_mda_y_2021 <- c(DP_july_21, DP_aug_21, DP_sept_21)
  DP_mda_y_2022 <- c(DP_july_22, DP_aug_22, DP_sept_22)
  
  
  #DP_mda_y_2022 <- c(y_2022+july, y_2022+aug, y_2022+sept)
  
  mda_events <- c(DP_mda_y_2021, DP_mda_y_2022)
  
  mda_distr <- length(mda_events)
  DP_cov <- unique(df_setting$mean_DP_cov)
  
  mdaparams_i <- set_mda(mdaparams_i,
                       drug = 1, #just DP
                       timesteps = mda_events,
                       coverages = rep(DP_cov, mda_distr), #modify cov for each distr acc to data
                       min_ages = rep(6*30, length(mda_events)), #min age is 6 months
                       max_ages = rep(85 * 365, length(mda_events))) #the max age of pop
  
  correlations <- get_correlation_parameters(mdaparams_i)
  correlations$inter_round_rho('bednets', 1) #same people get ITNs each year
  correlations$inter_intervention_rho('bednets', 'mda', -1) #not necessarily same people getting ITNs and DP
  
  # store in list
  sim_list[[i]] <- mdaparams_i
}

# assign scenario names after the loop
names(sim_list) <- paste0(
  "Q0_", df_scenario$Q0,
  "_phi_", df_scenario$bites_Bed,
  "_net_", df_scenario$retention_net_in_ 
)

#run the model: loop over sim_list and call run simulation 

run_control_list <- lapply(
  names(sim_list), 
  function(scenario){
    c("Running control scenario", scenario, "\n")
    run_simulation(
      timesteps = sim_length, 
      parameters = sim_list[[scenario]], 
      correlations = correlations
    )
  }
)

beepr::beep(sound = 1)

#names(run_control_list) <- paste0(
#  "Q0_", df_scenario$Q0,
#  "_phi_", df_scenario$bites_Bed,
#  "_net_", df_scenario$retention_net_in
#)
#ggplot(run_control_list$Q0_0.824_phi_0.836_net_360, aes(x = timestep, y = n_detect_pcr_0_29200/n_age_0_29200))+
#  geom_line()

# Assign numeric scenario names
names(run_control_list) <- paste0("Scenario_", seq_along(run_control_list))

# Check
names(run_control_list)
# [1] "Scenario_1" "Scenario_2" ... "Scenario_9"

results_df_control <- map2_dfr(
  run_control_list,
  names(run_control_list),
  ~ mutate(.x, scenario = .y)
)
ggplot(results_df_control, aes(x = ((timestep-1)/365)+2012, 
                       y = n_detect_pcr_0_29200 / n_age_0_29200, col = as.factor(scenario))) +
  geom_line() +
  #facet_wrap(~ scenario) +
  scale_x_continuous(breaks = 2012:2024) +
  coord_cartesian(xlim = c(2012, 2024 - 0.55)) +
  labs(x = "Year", y = "PCR prevalence")

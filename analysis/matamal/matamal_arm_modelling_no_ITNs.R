#how much does explicit modelling of ITNs influence impact estimates?
#repeat with an init_EIR of 7 but no ITNs this time. 


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

simparams1 <- get_parameters(overrides = list(
  human_population = human_population, 
  endec = TRUE, 
  bednets = TRUE,
  prevalence_rendering_min_ages = 0, 
  prevalence_rendering_max_ages = 80 * 365, #all age prevalence - in line with odin
  individual_mosquitoes = FALSE, 
  age_group_rendering_min_ages = age_min, 
  age_group_rendering_max_ages = age_max ,
  model_seasonality = TRUE,
  g0 = seasonality_params$g0,
  g = c(seasonality_params$g1, seasonality_params$g2, seasonality_params$g3),
  h = c(seasonality_params$h1, seasonality_params$h2, seasonality_params$h3)
  
))

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

demo_params <- add_demography(p = simparams1, demography=site_demo_2012)

mosq_params <- gamb_params 

#setting parameters
df_setting <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_odin_inputs.rds")

df_phi <- read.csv("analysis/matamal/phi_estimates_bijagos.csv", header = TRUE) %>%
  filter(scenario == "Data")


mosq_params$Q0 <- df_setting$hbi_median[1]
mosq_params$phi_bednets <- round(df_phi$mean,3) #phi-bed from Bijagos

simparams <- set_species(demo_params, species = list(mosq_params), 
                         proportions = c(1))

#set demography####


#set drugs and treatment. Use site file for AL####

#leaving out for now - can't find AL - it is case detection and treatment with AL and IPTp in pregnancy

#Update parameter set with chosen drug-specific parameters (AL and DHA/PQP)
drug_params <- set_drugs(simparams, list(AL_params))

sf_treatment <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/1.data-cleaning/output/site_file_tx_cov_AL.rds")
names(sf_treatment) <- c("year", "tx_cov")
sf_treatment_2012 <- sf_treatment %>%
  filter(year >= 2012)



# Set treatment program for AL (drug index = 1)
treatment_params <- set_clinical_treatment(
  parameters = drug_params,
  drug = 1, #just AL
  timesteps =  (sf_treatment_2012$year - 2012)*365 + 1, # Treatment coverage changes on day 300 and day 600
  coverages =  sf_treatment_2012$tx_cov) # The initial treatment coverage (0%) is the default. site file would multiply by prop_ACT and here we assume that to be 1  


##vector control####
#ITNs: setup


#retention_net <- c(12, 21, 38)*30 #AG work https://www.medrxiv.org/content/10.1101/2025.08.27.25334550v1.full.pdf
#
#df_nets <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_net_matamal_info.rds")
#
#df_nets2 <- df_nets %>%
#  filter(year >= 2012)
#
#itn_distr <- nrow(df_nets2)
#
#distr_campaign <- 6*30 # mass distributions are around June
#
#med_dn0 <- unique(df_nets2$dn0_med)
#med_rn0 <- unique(df_nets2$rn0_med)
#med_gamman <- unique(df_nets2$gamman_med)
#
##some sort of bug here in set_bednets
#bednet_params <- set_bednets(treatment_params, 
#                             timesteps = (df_nets2$year - 2012)*365 + distr_campaign, 
#                             retention = retention_net[2], 
#                             coverages = df_nets2$itn_input_distr,
#                             dn0 = matrix(rep(med_dn0, itn_distr), nrow = itn_distr,ncol = 1), # Matrix of death probabilities for each mosquito species over time
#                             rn = matrix(rep(med_rn0, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
#                             rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
#                             gamman = rep(med_gamman * 365, itn_distr))
#
#correlations <- get_correlation_parameters(bednet_params)
#correlations$inter_round_rho('bednets', 1)


#use cali or malariaeqm methods to get the correct starting_EIR
#how do you specify that this prev survey is from 2019


#nov <- 30*11
#y_2018 <- 365*7 #7 into simulation
#prev_survey_date <- nov+y_2018


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


prev_survey_date <- timestep_from_2012(2018, 11, 1) 



baseline_prev_2018 <- 0.148 # from DTNmapper
################################################################################
library(cali)

# Prepare a summary function that returns the mean PfPR2-10 from each simulation output: 
summary_mean_pfpr_all_age <- function (x) {
  
  # Calculate the PfPR2-10:
  prev <- x$n_detect_pcr_0_29200[prev_survey_date]/x$n_age_0_29200[prev_survey_date] #specify here that this is around Nov 2018
  
  # Return the calculated PfPR all age:
  return(prev)
}

# Establish a target PfPR2-10 value:
target_pfpr <- baseline_prev_2018

treatment_params$timesteps <- sim_length


# Run the calibrate() function:
cali_EIR <- cali::calibrate(target = target_pfpr,
                            eq_prevalence = 0.1, #just a starting point for search
                            summary_function = summary_mean_pfpr_all_age,
                            human_population = c(1e+04, 1e+05, 1e+06),
                            parameters = treatment_params, 
                            eir_limits = c(2,12)) 
#cali_EIR_in <- 6.38
#check cali output#########################
#
##try some trial runs with a slightly lower EIR
#
###############################################################################################
#bednet_params$timesteps <- sim_length
##changing caliEIR, currently 9.11..
#simparams_cali <- set_equilibrium(bednet_params, init_EIR = 7)
#
## Run the simulation:
#cali_sim <- run_simulation(timesteps = (simparams_cali$timesteps),
#                           parameters = simparams_cali, 
#                           correlations = correlations)
#
#ggplot(cali_sim, aes(x = (timestep/365), y = n_detect_pcr_0_29200/n_age_0_29200))+
#  geom_line()+
#  geom_vline(aes(xintercept = (prev_survey_date/365)), col = "red")
#
#
#ggplot(cali_sim, aes(x = ((timestep+1)/365)+2012, y = n_detect_pcr_0_29200/n_age_0_29200))+
#  geom_line()+
#  geom_vline(aes(xintercept = (prev_survey_date/365)+2012), col = "red")+
#  scale_x_continuous(breaks = 2012:2025)
#
#cali_sim$n_detect_pcr_0_29200[prev_survey_date]/cali_sim$n_age_0_29200[prev_survey_date] #15.04%
#cali_sim$n_detect_lm_0_29200[prev_survey_date]/cali_sim$n_age_0_29200[prev_survey_date] #~0.07801603 all age prev



input_EIR <- 7
#############################################
#then set eqm
simparams <- set_equilibrium(parameters = treatment_params, 
                             init_EIR = input_EIR)



#then put in the DP MDA
# Make a copy of the base simulation parameters to which we can add the MDA campaign parameters:
mdaparams <- simparams

# Update the parameter list with the default parameters for sulphadoxine-pyrimethamine
# amodiaquine (SP-AQ)
mdaparams <- set_drugs(mdaparams, list(DHA_PQP_params))

# Specify the days on which to administer: modify according to the trial information
# for purpose here, could do first of Jul, Aug and Sept of 2021 and 2022

#july <- 30*7
#aug <- 30*8
#sept <- 30*9
#
#y_2021 <- 365*10
#y_2022 <- 365*11 

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

# Use set_mda() function to set the proportion of the population that the MDA reaches and the age 
# ranges which can receive treatment.
mdaparams <- set_mda(mdaparams,
                     drug = 1, #just DP
                     timesteps = mda_events,
                     coverages = rep(DP_cov, mda_distr), #modify cov for each distr acc to data
                     min_ages = rep(6*30, length(mda_events)), #min age is 6 months
                     max_ages = rep(85 * 365, length(mda_events))) #the max age of pop

##################################################################################################
#run simulation here for control arm
run_control <- run_simulation(timesteps = sim_length,
                              parameters = mdaparams) 

run_control <- run_control %>%
  mutate(year = ((timestep - 1) / 365) + 2012, 
         year_sim = ((timestep - 1) / 365))

run_control %>%
  filter(timestep == 730) %>%
  select(timestep, year)

run_control 

#prevalence estimates from trial to overlay 
prev<- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/baseline_prev_malariasim_params.rds")
prev <- prev %>%
  add_row(year = as.factor(2018), mean_prev = 0.148, lower_prev = 14.5, upper_prev = 14.9)



#DP_mda_y_2021 <- c(y_2021+july, y_2021+aug, y_2021+sept)
#
#DP_mda_y_2022 <- c(y_2022+july, y_2022+aug, y_2022+sept)

#prev survey 2021 31st Oct - 22nd Nov 2021. Just say 1st Nov
#prev survey 2022 4th Nov - 1st Dec 2022
nov <- 30*11
#y_2019 <- 365*8 # 8 years into simulation
#prev_survey_date <- nov+y_2018

#y_2021 <- 365*10
prev_survey_date_2021 <- timestep_from_2012(2021, 11, 1)

#y_2022 <- 365*11
prev_survey_date_2022 <- timestep_from_2012(2022, 11, 1)

prev_2018 <- prev %>%
  filter(year == 2018) %>%
  mutate(year_sim = 6,
         day = prev_survey_date) #for 2018

prev_2021 <- prev %>%
  filter(year == 2021) %>%
  mutate(day = prev_survey_date_2021, 
         year_sim = 9)

prev_2022 <- prev %>%
  filter(year == 2022) %>%
  mutate(day = prev_survey_date_2022, 
         year_sim = 10)



#now put this on year scale

#add the DP arrows

df_prev_cluster <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/cluster_prev_clean.rds")

df_prev_cluster <- df_prev_cluster %>%
  mutate(day = case_when(year == 2021 ~ prev_2021$day, 
                         year == 2022 ~ prev_2022$day, 
                         TRUE ~ NA_real_))

df_prev_cluster_control <- df_prev_cluster %>%
  filter(arm == 0)

df_prev_cluster_int <- df_prev_cluster %>%
  filter(arm == 1)


distr_campaign <- 6*30 #distribution in June
distr_years <- c(2014, 2017, 2020)



ggplot(run_control, aes(x = (timestep/365), y = (n_detect_pcr_0_29200/n_age_0_29200)))+
  geom_line()+
  scale_x_continuous(breaks = 0:14)
######################

mda_campaign_arrows <- data.frame(
  x = ((mda_events-1)/365)+2012,
  xend = ((mda_events-1)/365)+2012,
  y = rep(0.2, length(distr_years)),
  yend = rep(0.17, length(distr_years))
)

#saveRDS(run_control, "analysis/matamal/output/run_control_DP.rds")

control_plot <- ggplot(run_control, aes(x = ((timestep-1)/365)+2012, y = (n_detect_pcr_0_29200/n_age_0_29200)))+
  geom_line()+
  #coord_cartesian(xlim = c(2017, 2023))+
  #ITN campaigns 
  annotate("segment", x = distr_years[1]+(distr_campaign/365), 
           xend = distr_years[1]+(distr_campaign/365), y = 0.2, yend = 0.17, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[2]+(distr_campaign/365), 
           xend = distr_years[2]+(distr_campaign/365), y = 0.2, yend = 0.17, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[3]+(distr_campaign/365), 
           xend = distr_years[3]+(distr_campaign/365), y = 0.2, yend = 0.17,
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  #2018 survey
  geom_point(data = prev_2018, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "red", shape = 4, 
             inherit.aes = FALSE)+
  geom_errorbar(data = prev_2018,
                aes(x = ((day-1)/365)+2012, ymin = lower_prev/100, ymax = upper_prev/100),
                color = "red", width = 0.2,
                inherit.aes = FALSE)+
  coord_cartesian(xlim = c(2012, 2023 -0.5)) +
  scale_x_continuous(breaks = 2012:2022)+ #no ITN inputs in df for 2023
  
  #2021 survey
  geom_point(data = prev_2021, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "red", shape = 4, 
             inherit.aes = FALSE)+
  geom_errorbar(data = prev_2021,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE)+
  #2022 survey
  geom_errorbar(data = prev_2022,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE)+
  geom_point(data = prev_2022, aes(x = (day/365)+2012, y = mean_prev), col = "red", shape = 4, 
             inherit.aes = FALSE)+
  
  #cluster level data 
  geom_point(data = df_prev_cluster_control, aes(x = (day/365)+2012, y = all_age_prev/100), col = "slateblue")+
  
  theme_bw()+
  labs(x = "Year", y = "PCR prevalence")+
  ggtitle("Control arm")+
  geom_segment(
    data = mda_campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1.5, col = "orange"
  )

campaign_arrows <- data.frame(
  x = distr_years + (distr_campaign / 365),
  xend = distr_years + (distr_campaign / 365),
  y = rep(0.80, length(distr_years)),
  yend = rep(0.75, length(distr_years))
)

ggplot(run_control, aes(x = (timestep/365) + 2012, y = (n_use_net / human_population))) +
  geom_line() +
  coord_cartesian(ylim = c(0, 1)) +
  geom_segment(
    data = campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1.5
  ) +
  labs(x = "Year", y = "Proportion using net") +
  theme_bw()+
  scale_x_continuous(breaks = 2012:2022)



#for intervention arm, add endectocides in 

#endectocide set up####
mda_int <- 30
eff_len <- 23
#early_IVM <- 180
#endec_y1_start <- july + y_2021 #start time of first MDA
#endec_y2_start <- july+ y_2022

endec_y1_start <- timestep_from_2012(2021, 7, 1)
endec_y2_start <- timestep_from_2012(2022, 7, 1)


#start time of each MDA (Jul, Aug, Sep 2021 and Jul, Aug, Sep 2022)
early_IVM_start <- c(endec_y1_start, endec_y1_start+mda_int, endec_y1_start+mda_int+mda_int, 
                     endec_y2_start, endec_y2_start+mda_int, endec_y2_start+mda_int+mda_int)
#overall time that endec is "on"
early_steps <- (early_IVM_start[1]):(early_IVM_start[6]+eff_len)

#actual on period for each MDA..so increase this up to 6 (3 per year)
early_endec_on <- c(as.numeric(early_IVM_start[1]):as.numeric(early_IVM_start[1]+eff_len), 
                    as.numeric(early_IVM_start[2]):as.numeric(early_IVM_start[2]+eff_len), 
                    as.numeric(early_IVM_start[3]):as.numeric(early_IVM_start[3]+eff_len), 
                    as.numeric(early_IVM_start[4]):as.numeric(early_IVM_start[4]+eff_len),
                    as.numeric(early_IVM_start[5]):as.numeric(early_IVM_start[5]+eff_len),
                    as.numeric(early_IVM_start[6]):as.numeric(early_IVM_start[6]+eff_len)
)


#start time of each MDA for the model input
early_endec_ts <- c(early_IVM_start[1], early_IVM_start[2], early_IVM_start[3], 
                    early_IVM_start[4], early_IVM_start[5], early_IVM_start[6])



#early MDA with long net retention
endec_params <- set_endectocide(mdaparams, timesteps = early_steps,
                                endec_on = early_endec_on, 
                                endec_ts = early_endec_ts)

#run model####
run_endec <- run_simulation(endec_params, 
                            timesteps = sim_length)


#run_simulation(endec_params,  timesteps = sim_length,correlations = correlations_long)

#saveRDS(run_endec, file = "analysis/matamal/output/run_int_DP_endec.rds")

int_plot <- ggplot(run_endec, aes(x = ((timestep-1)/365)+2012, y = (n_detect_pcr_0_29200/n_age_0_29200)))+
  geom_line()+
  
  #2019 survey
  geom_point(data = prev_2018, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "red", shape = 4, 
             inherit.aes = FALSE)+
  geom_errorbar(data = prev_2018,
                aes(x = ((day-1)/365)+2012, ymin = lower_prev/100, ymax = upper_prev/100),
                color = "red", width = 0.2,
                inherit.aes = FALSE)+
  #2021 survey
  geom_point(data = prev_2021, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "red", shape = 4, 
             inherit.aes = FALSE)+
  geom_errorbar(data = prev_2021,
                aes(x = ((day-1)/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE)+
  #2022 survey
  geom_errorbar(data = prev_2022,
                aes(x = ((day-1)/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE)+
  geom_point(data = prev_2022, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "red", shape = 4, 
             inherit.aes = FALSE)+
  #cluster level data 
  geom_point(data = df_prev_cluster_int, aes(x = ((day-1)/365)+2012, y = all_age_prev/100), col = "slateblue")+
  coord_cartesian(xlim = c(2012, 2023 -0.5)) +
  scale_x_continuous(breaks = 2012:2022)   + #no ITN inputs in df for 2023
  theme_bw()+
  labs(x = "Year", y = "PCR prevalence")+
  ggtitle("Intervention arm")+
  #scale_x_continuous(breaks = 2012:2023)+
  #ITN campaigns 
  annotate("segment", x = distr_years[1]+(distr_campaign/365), 
           xend = distr_years[1]+(distr_campaign/365), y = 0.2, yend = 0.17, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[2]+(distr_campaign/365), 
           xend = distr_years[2]+(distr_campaign/365), y = 0.2, yend = 0.17, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[3]+(distr_campaign/365), 
           xend = distr_years[3]+(distr_campaign/365), y = 0.2, yend = 0.17,
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  geom_segment(
    data = mda_campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1.5, col = "orange"
  )

cowplot::plot_grid(control_plot, int_plot)

#actually do repeat the cali step - get to the baseline survey prevalence without ITNs explicitly modelled- what Hannah did

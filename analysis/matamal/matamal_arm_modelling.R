devtools::load_all()
require(tidyverse)
require(cali)
#malariasim with baseline scenario, ITN and endectocides

#baseline scenario####

#matamal trial modelling

year <- 365
sim_length <- 365*11 # 10 years. From 2012 to 2022
human_population <- 1e5
#starting_EIR <- 10


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
  age_group_rendering_max_ages = age_max ,
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

#set demography####


#set drugs and treatment. Use site file for AL####

#leaving out for now - can't find AL - it is case detection and treatment with AL and IPTp in pregnancy

# Update parameter set with chosen drug-specific parameters (AL and DHA/PQP)
#drug_params <- set_drugs(simparams, list(AL_params))
#
#
## Set treatment program for AL (drug index = 1)
#treatment_params <- set_clinical_treatment(
#  parameters = drug_params,
#  drug = 1, #just AL
#  timesteps =  c(300,600), # Treatment coverage changes on day 300 and day 600
#  coverages =  c(0.4,0)) # The initial treatment coverage (0%) is the default 
## and does not need to be set
#
#
##simparams <- set_equilibrium(parameters = treatment_params, 
# #                            init_EIR = starting_EIR)

##vector control####
#ITNs: setup


retention_net <- 21*30 # 730 days. Andrew Glover work.

df_nets <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_net_matamal_info.rds")

df_nets2 <- df_nets %>%
  filter(year >= 2012)

itn_distr <- nrow(df_nets2)

distr_campaign <- 6*30 # mass distributions are around June

med_dn0 <- unique(df_nets2$dn0_med)
med_rn0 <- unique(df_nets2$rn0_med)
med_gamman <- unique(df_nets2$gamman_med)

#some sort of bug here in set_bednets
bednet_params <- set_bednets(simparams, 
                             timesteps = (df_nets2$year - 2012)*365 + distr_campaign, 
                             retention = retention_net, 
                             coverages = df_nets2$itn_input_distr,
                             dn0 = matrix(rep(med_dn0, itn_distr), nrow = itn_distr,ncol = 1), # Matrix of death probabilities for each mosquito species over time
                             rn = matrix(rep(med_rn0, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
                             rnm = matrix(rep(0.24, itn_distr), nrow = itn_distr, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
                             gamman = rep(med_gamman * 365, itn_distr))

correlations <- get_correlation_parameters(bednet_params)
correlations$inter_round_rho('bednets', 1)


#use cali or malariaeqm methods to get the correct starting_EIR
#how do you specify that this prev survey is from 2019


nov <- 30*11
y_2019 <- 365*7 # 7 years into simulation
prev_survey_date <- nov+y_2019

baseline_prev_2019 <- unique(df_setting$mean_prev) #from late October to early December 2019

library(cali)

# Prepare a summary function that returns the mean PfPR2-10 from each simulation output: 
summary_mean_pfpr_all_age <- function (x) {
  
  # Calculate the PfPR2-10:
  prev <- x$n_detect_pcr_0_31025[prev_survey_date]/x$n_age_0_31025[prev_survey_date] #specify here that this is around Nov 2019
  
  # Return the calculated PfPR2-10:
  return(prev)
}

# Establish a target PfPR2-10 value:
target_pfpr <- baseline_prev_2019

# Add a parameter to the parameter list specifying the number of timesteps to
# simulate over. Note, increasing the number of steps gives the simulation longer
# to stablise/equilibrate, but will increase the runtime for calibrate().  
#simparams$timesteps <- 3 * 365

# Establish a tolerance value:
pfpr_tolerance <- 0.01

# Set upper and lower EIR bounds for the calibrate function to check (remembering EIR is
# the variable that is used to tune to the target PfPR):
lower_EIR <- 5; upper_EIR <- 15

# Run the calibrate() function:
cali_EIR <- calibrate(target = target_pfpr,
                      summary_function = summary_mean_pfpr_all_age,
                      parameters = simparams,
                      tolerance = pfpr_tolerance, 
                      low = lower_EIR, high = upper_EIR) #adjust this

simparams_cali <- set_equilibrium(simparams, init_EIR = cali_EIR)

# Run the simulation:
cali_sim <- run_simulation(timesteps = (simparams_cali$timesteps),
                           parameters = simparams_cali)



cali_pfpr <- cali_sim$n_detect_pcr_0_31025 / cali_sim$n_age_0_31025 


# Set the plotting window:
par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))

# Plot the PfPR2-10 under the EIR recommended by the cali method:
plot(x = df$timestep,
     y = df$malsim_pfpr,
     type = "b",
     ylab = expression(paste(italic(Pf),"PR"[2-10])),
     xlab = "Time (days)",
     ylim = c(target_pfpr - 0.2, target_pfpr + 0.2),
     #ylim = c(0, 1),
     col = cols[3])

# Add grid lines
grid(lty = 2, col = "grey80", lwd = 0.5)

# Add a textual identifier and lines indicating target PfPR2-10 with
# tolerance bounds:
text(x = 10, y = 0.47, pos = 4, cex = 0.9,
     paste0("a) cali \n init_EIR = ", round(cali_EIR, digits = 3)))
abline(h = target_pfpr, col = "dodgerblue", lwd = 2)
abline(h = target_pfpr - pfpr_tolerance, lwd = 2)
abline(h = target_pfpr + pfpr_tolerance, lwd = 2)

# Plot the PfPR2-10 under the EIR recommended by the malariasimulation method:
plot(x = df$timestep,
     y = df$cali_pfpr,
     type = "b",
     xlab = "Time (days)",
     ylab = "",
     ylim = c(target_pfpr - 0.2, target_pfpr + 0.2),
     #ylim = c(0, 1),
     col = cols[1])

# Add grid lines
grid(lty = 2, col = "grey80", lwd = 0.5)

# Add a textual identifier and lines indicating target PfPR2-10 with
# tolerance bounds
text(x = 10, y = 0.47, pos = 4, cex = 0.9,
     paste0("b) malariasimulation \n init_EIR = ", round(malsim_EIR, digits = 3)))
abline(h = target_pfpr, col = "dodgerblue", lwd = 2)
abline(h = target_pfpr - pfpr_tolerance, lwd = 2)
abline(h = target_pfpr + pfpr_tolerance, lwd = 2)






#then set eqm
simparams <- set_equilibrium(parameters = simparams, 
                             init_EIR = starting_EIR)

out <- run_simulation(simparams, timesteps = sim_length)

ggplot(out, aes(x = timestep, y = n_detect_lm_0_31025/n_age_0_31025))+
  geom_line()+
  ylim(0,1)

#
#then put in the DP MDA
# Make a copy of the base simulation parameters to which we can add the MDA campaign parameters:
mdaparams <- simparams

# Update the parameter list with the default parameters for sulphadoxine-pyrimethamine
# amodiaquine (SP-AQ)
mdaparams <- set_drugs(mdaparams, list(DHA_PQP_params))

# Specify the days on which to administer: modify according to the trial information
# for purpose here, could do first of Jul, Aug and Sept
mda_events <- (c(1, 2) * 365)

# Use set_mda() function to set the proportion of the population that the MDA reaches and the age 
# ranges which can receive treatment.
mdaparams <- set_mda(mdaparams,
                     drug = 1, #just DP
                     timesteps = mda_events,
                     coverages = rep(.8, 2), #modify cov for each distr acc to data
                     min_ages = rep(0, length(mda_events)), #check minimum age for DP
                     max_ages = rep(200 * 365, length(mda_events)))


#run simulation here for control arm

#for intervention arm, add endectocides in 

#endectocide set up####
mda_int <- 30
eff_len <- 23
#early_IVM <- 180
endec_y1_start <- 365 #start time of first MDA
endec_y2_start <- 365*2
#start time of each MDA (Jul, Aug, Sep 2021 and Jul, Aug, Sep 2022)
early_IVM_start <- c(endec_y1_start, endec_y1_start+mda_int, endec_y1_start+mda_int+mda_int, 
                     endec_y2_start, endec_y2_start+mda_int, endec_y2_start+mda_int+mda_int)

#start time of each MDA for the model input
early_endec_ts <- c(early_IVM_start[1], early_IVM_start[2], early_IVM_start[3], 
                    early_IVM_start[4], early_IVM_start[5], early_IVM_start[6])

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

#early MDA with long net retention
endec_params <- set_endectocide(simparams, timesteps = early_steps,
                                endec_on = early_endec_on, 
                                endec_ts = early_endec_ts)

#run model####
run_endec <- run_simulation(endec_params, 
                            timesteps = sim_length)
                            #correlations = correlations_long)
                            
#run_simulation(endec_params,  timesteps = sim_length,correlations = correlations_long)

ggplot(run_endec, aes(x = timestep, y = (n_detect_lm_0_31025/n_age_0_31025)))+
  geom_line()+
  geom_vline(aes(xintercept = early_IVM_start[3]+30+30), lty = "dashed")+
  geom_vline(aes(xintercept = early_IVM_start[6]+30+30), lty = "dashed")


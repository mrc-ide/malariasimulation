# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# #+++++ antimalarial_resistance/etf development script +++++#
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# 
# ##' README:
# ##' THis script is designed as a place to develop and store the changes being implemented into the
# ##' malariasimulation package in support of integrating antimalarial resistance. This file is
# ##' designed to be deleted prior to any merge and is purely for development purposes.
# 
# ##' CURENTLY WORKING ON: amending the calculate_treated() function to update
# ##'
# ##'
# ##' Use Ctrl + Shft + C to comment everything in/out for load_all() to work
# 
# #----- 1) Preamble ---------------------------------------------------------------------------------
# 
# # Load in the individual package:
# library(individual)
# library(devtools)
# 
# # Load the malariasimualtion functions (make sure this script is commented out when you run this)
# devtools::load_all(".")
# 
# #----- 2) Original calculate_treated() function ----------------------------------------------------
# 
# ##' This is the original calculate_treated() function, renamed to avoid any conflicts.
# calculate_treated <- function(
#     variables,
#     clinical_infections,
#     parameters,
#     timestep,
#     renderer
# ) {
# 
#   # Gather the treatment coverage(s) in the current timestep:
#   treatment_coverages <- get_treatment_coverages(parameters, timestep)
# 
#   # Sum to get the total treatment coverage in the current time step:
#   ft <- sum(treatment_coverages)
# 
#   # If treatment = 0, return an empty Bitset for the treated individuals:
#   if (ft == 0) {
#     return(individual::Bitset$new(parameters$human_population))
#   }
# 
#   # Add the total treatment coverage to the renderer:
#   renderer$render('ft', ft, timestep)
# 
#   # Sample individuals from clinically infected to seek treatment using treatment coverage:
#   seek_treatment <- sample_bitset(clinical_infections, ft)
# 
#   # Store the number of people that seek treatment this time step:
#   n_treat <- seek_treatment$size()
# 
#   # Add the number of people seeking treatment in the current time step to the renderer:
#   renderer$render('n_treated', n_treat, timestep)
# 
#   # Assign each individual seeking treatment a drug based on their coverage(s):
#   drugs <- as.numeric(parameters$clinical_treatment_drugs[
#     sample.int(
#       length(parameters$clinical_treatment_drugs),
#       n_treat,
#       prob = treatment_coverages,
#       replace = TRUE
#     )
#   ])
# 
#   # Determine, using drug efficacies, use bernoulli trials to see who's treatment worked:
#   successful <- bernoulli_multi_p(parameters$drug_efficacy[drugs])
# 
#   # Subset the successfully treated people:
#   treated_index <- bitset_at(seek_treatment, successful)
# 
#   # For those people who have been successfully treated:
#   if (treated_index$size() > 0) {
# 
#     # Queue update to the infectious state to Tr:
#     variables$state$queue_update('Tr', treated_index)
# 
#     # Queue update to the infectivity to reflect their new state (Tr):
#     variables$infectivity$queue_update(
#       parameters$cd * parameters$drug_rel_c[drugs[successful]],
#       treated_index
#     )
#     # Queue update to the last drug each successfully treated person received:
#     variables$drug$queue_update(
#       drugs[successful],
#       treated_index
#     )
#     # Queue update to the time step on which each successfully treated person last
#     # received a drug:
#     variables$drug_time$queue_update(
#       timestep,
#       treated_index
#     )
#   }
# 
#   # Return the Bitset of treated individuals:
#   treated_index
# }
# 
# #----- 3) Function Development Set-up --------------------------------------------------------------
# 
# # Set an base parameter list
# simparams <- get_parameters()
# 
# # Add AL and DHA_PQP drug parameters to parameter list
# simparams <- set_drugs(simparams, list(AL_params, DHA_PQP_params))
# 
# # Add clinical treatment events using set_clinical_treatment():
# simparams <- set_clinical_treatment(parameters = simparams, drug = 2,
#                                     timesteps = 25, coverages = 0.2)
# simparams <- set_clinical_treatment(parameters = simparams, drug = 1,
#                                     timesteps = 50, coverages = 0.4)
# 
# # Check the clinical treatment order and coverages:
# simparams$clinical_treatment_drugs
# simparams$clinical_treatment_coverages
# 
# ##' AL is drug 1, used second clinically at a coverage of 0.4
# ##' DHA-PQP is drug 2, used first clinically at a coverage of 0.2
# 
# # Add resistance using set_antimalarial_resistance() using same order as clinical treatment
# simparams <- set_antimalarial_resistance(parameters = simparams,
#                                          drug = 2,
#                                          timesteps = c(0, 25),
#                                          artemisinin_resistance = c(0.3, 0.42),
#                                          partner_drug_resistance = c(0, 0),
#                                          slow_parasite_clearance_prob = rep(0.1, 2),
#                                          early_treatment_failure_prob = rep(0.6, 2),
#                                          late_clinical_failure_prob = rep(0.1, 2),
#                                          late_parasitological_prob = rep(0.1, 2),
#                                          reinfection_prob = rep(0.1, 2))
# simparams <- set_antimalarial_resistance(parameters = simparams,
#                                          drug = 1,
#                                          timesteps = c(0, 50),
#                                          artemisinin_resistance = c(0, 0),
#                                          partner_drug_resistance = c(0, 0.5),
#                                          slow_parasite_clearance_prob = rep(0.05, 2),
#                                          early_treatment_failure_prob = rep(0.5, 2),
#                                          late_clinical_failure_prob = rep(0.5, 2),
#                                          late_parasitological_prob = rep(0.5, 2),
#                                          reinfection_prob = rep(0.5, 2))
# 
# # Check the resistance parameters:
# simparams$antimalarial_resistance_drug
# simparams$antimalarial_resistance_timesteps
# simparams$prop_artemisinin_resistant; simparams$prop_partner_drug_resistant
# simparams$early_treatment_failure_prob; simparams$slow_parasite_clearance_prob
# simparams$late_clinical_failure_prob; simparams$late_parasitological_failure_prob; simparams$reinfection_during_prophylaxis
# 
# # Set up the model variables:
# variables <- create_variables(parameters = simparams)
# 
# # Establish Bitset of the human population
# current_pop <- Bitset$new(simparams$human_population, from = NULL)
# current_pop$insert(1:simparams$human_population); current_pop$to_vector()
# 
# # Subset a group of clinically infected individuals from the population
# clinical_infections <- sample_bitset(b = current_pop, rate = 0.4)
# clinical_infections$to_vector(); clinical_infections$size()
# 
# ##' TIMESTEP
# timesteps <- 100
# current_timestep <- 50
# 
# ##' RENDERER
# renderer <- individual::Render$new(timesteps)
# renderer$to_dataframe()
# 
# #----- 4) Individual Steps (Original) --------------------------------------------------------------
# 
# # Gather the treatment coverage(s) in the current timestep:
# treatment_coverages <- get_treatment_coverages(simparams, current_timestep)
# 
# # Sum to get the total treatment coverage in the current time step:
# ft <- sum(treatment_coverages)
# 
# # If treatment = 0, return an empty Bitset for the treated individuals:
# if (ft == 0) {
#   return(individual::Bitset$new(simparams$human_population))
# } else {
#   print("OK")
# }
# 
# # Add the total treatment coverage to the renderer:
# renderer$render('ft', ft, current_timestep)
# renderer$to_dataframe()[current_timestep,]
# 
# # Sample individuals from clinically infected to seek treatment using treatment coverage:
# seek_treatment <- sample_bitset(clinical_infections, ft)
# seek_treatment$to_vector(); seek_treatment$size()
# 
# # Store the number of people that seek treatment this time step:
# n_treat <- seek_treatment$size(); n_treat
# 
# # Add the number of people seeking treatment in the current time step to the renderer:
# renderer$render('n_treated', n_treat, current_timestep)
# renderer$to_dataframe()[current_timestep,]
# 
# # Assign each individual seeking treatment a drug based on their coverage(s):
# drugs <- as.numeric(simparams$clinical_treatment_drugs[
#   sample.int(
#     length(simparams$clinical_treatment_drugs),
#     n_treat,
#     prob = treatment_coverages,
#     replace = TRUE
#   )
# ])
# drugs; length(drugs) == n_treat
# 
# # Determine, using drug efficacies, use bernoulli trials to see who's treatment worked:
# successful <- bernoulli_multi_p(simparams$drug_efficacy[drugs])
# successful; length(successful); (length(successful)/length(drugs))*100
# 
# # Subset the successfully treated people:
# treated_index <- bitset_at(seek_treatment, successful)
# treated_index$to_vector(); treated_index$size()
# 
# 
# # # For those people who have been successfully treated:
# # if (treated_index$size() > 0) {
# #
# #   # Queue update to the infectious state to Tr:
# #   variables$state$queue_update('Tr', treated_index)
# #
# #   # Queue update to the infectivity to reflect their new state (Tr):
# #   variables$infectivity$queue_update(
# #     parameters$cd * parameters$drug_rel_c[drugs[successful]],
# #     treated_index
# #   )
# #   # Queue update to the last drug each successfully treated person received:
# #   variables$drug$queue_update(
# #     drugs[successful],
# #     treated_index
# #   )
# #   # Queue update to the time step on which each successfully treated person last
# #   # received a drug:
# #   variables$drug_time$queue_update(
# #     timestep,
# #     treated_index
# #   )
# # }
# 
# # Return the Bitset of treated individuals:
# treated_index
# 
# #----- 5) Individual Steps (Updated) ---------------------------------------------------------------
# 
# # Gather the treatment coverage(s) in the current timestep:
# treatment_coverages <- get_treatment_coverages(simparams, current_timestep); treatment_coverages
# 
# # Sum to get the total treatment coverage in the current time step:
# ft <- sum(treatment_coverages); ft
# 
# # If treatment = 0, return an empty Bitset for the treated individuals:
# if (ft == 0) {
#   return(individual::Bitset$new(simparams$human_population))
# } else {
#   print("OK")
# }
# 
# # Add the total treatment coverage to the renderer:
# renderer$render('ft', ft, current_timestep)
# renderer$to_dataframe()[current_timestep,]
# 
# # Sample individuals from clinically infected to seek treatment using treatment coverage:
# seek_treatment <- sample_bitset(clinical_infections, ft)
# seek_treatment$to_vector(); seek_treatment$size()
# paste0(((seek_treatment$size()/clinical_infections$size())*100),
#        "% of clinically infected individuals sought treatment")
# 
# # Store the number of people that seek treatment this time step:
# n_treat <- seek_treatment$size(); n_treat
# 
# # Add the number of people seeking treatment in the current time step to the renderer:
# renderer$render('n_treated', n_treat, current_timestep)
# renderer$to_dataframe()[current_timestep,]
# 
# # Assign each individual seeking treatment a drug based on their coverage(s):
# drugs <- as.numeric(simparams$clinical_treatment_drugs[
#   sample.int(
#     length(simparams$clinical_treatment_drugs),
#     n_treat,
#     prob = treatment_coverages,
#     replace = TRUE
#   )
# ]); drugs; length(drugs) == n_treat; table(drugs)
# 
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# #+++++ RESISTANCE
# 
# # Set the resistance drug index (each individuals drug):
# AM_drug_index <- as.numeric(simparams$antimalarial_resistance_drug[drugs])
# drugs; AM_drug_index
# 
# # Generate a vector containing the proportion of artemisinin resistance to the drug taken by each
# # individual and it's probability of manifesting as early treatment failure:
# art_resistance <- vector(); etf_prob <- vector()
# for(i in 1:n_treat) {
# 
#   # Identify the correct resistance index:
#   last_resistance <- match_timestep(ts = simparams$antimalarial_resistance_timesteps[[AM_drug_index[i]]],
#                                     t = current_timestep)
# 
#   # Retrieve the resistance value corresponding to the drug and time step:
#   art_resistance[i] <- simparams$prop_artemisinin_resistant[[AM_drug_index[i]]][last_resistance]
# 
#   # Retrieve the ETF phenotype probability corresponding to the drug and time step:
#   etf_prob[i] <- simparams$early_treatment_failure_prob[[AM_drug_index[i]]][last_resistance]
# 
# }
# drugs; art_resistance; etf_prob
# 
# # The probability that an individual fails is the product of the proportion of individuals with resistance
# # to the drug they've taken and the probability that early treatment failure manifests in cases of resistance
# # to that drug. The probability individuals are treated successfully is therefore 1 - this product.
# 
# # Run bernoulli trials to determine which individuals fail treatment due to resistance to the drug
# # they have taken:
# susceptible <- bernoulli_multi_p(p = 1 - (art_resistance * etf_prob))
# length(susceptible); (length(susceptible)/n_treat)*100
# 
# ##' Drug 1 is AL and Drug 2 is DHA-PQP. Drug 2 is set first in clinical treatment and resistance.
# ##' The proportion of individuals in the current timestep with artemisinin resistance to drug 1 is
# ##' 0, while the proportion of the population with artemisinin resistance to drug 2 is 0.42. The prob.
# ##' of early treatment failure for Drug 1 is 0.5, and is 0.6 for Drug 2.
# 
# # Calculate the number of individuals who failed treatment due to early treatment failure:
# n_ETF <- n_treat - length(susceptible)
# 
# # Add the number of Early Treatment Failures to the renderer:
# renderer$render('n_ETF', n_ETF, current_timestep)
# renderer$to_dataframe()[current_timestep,]
# 
# # Remove the individuals who failed treatment due to resistance from the drugs vector:
# drugs[susceptible]; length(drugs[susceptible]); (length(drugs[susceptible])/n_treat)*100
# drugs_2 <- drugs[susceptible] 
# 
# # Check which individuals failed and which drug they took (only modelled resistance to drug 2) 
# res_fail_individuals <- which(!((1:n_treat) %in% susceptible))
# drugs[res_fail_individuals]; table(drugs); table(drugs_2)
# drugs; drugs_2
# 
# # Subset the people who remained susceptible
# seek_treatment$to_vector(); seek_treatment$size()
# susceptible_index <- bitset_at(seek_treatment, susceptible)
# susceptible_index$to_vector(); susceptible_index$size()
# 
# # So, we should now have a slightly shorter vector of drugs taken by each individual to account for
# # the individuals who gave failed treatment due to resistance. We can then use this to work out who
# # is successfully treated given the drug efficacy.
# 
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# #+++++ DRUG EFFICACY
# # Determine, using drug efficacies, use bernoulli trials to see who's treatment worked:
# successful <- bernoulli_multi_p(simparams$drug_efficacy[drugs_2])
# successful; length(successful); (length(successful)/length(drugs_2))*100
# 
# # Calculate the number of people who failed treatment due to efficacy:
# n_treat_eff_fail <- length(susceptible) - length(successful)
# 
# # Add the number of treatment efficacy failures to the renderer:
# renderer$render('n_treat_eff_fail', n_treat_eff_fail, current_timestep)
# renderer$to_dataframe()[current_timestep,]
# 
# # Subset the successfully treated people:
# treated_index <- bitset_at(susceptible_index, successful)
# treated_index$to_vector(); treated_index$size()
# 
# # Compare the three bitsets:
# seek_treatment$to_vector(); seek_treatment$size()
# susceptible_index$to_vector(); susceptible_index$size()
# treated_index$to_vector(); treated_index$size()
# 
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# 
# # # For those people who have been successfully treated:
# # if (treated_index$size() > 0) {
# #
# #   # Queue update to the infectious state to Tr:
# #   variables$state$queue_update('Tr', treated_index)
# #
# #   # Queue update to the infectivity to reflect their new state (Tr):
# #   variables$infectivity$queue_update(
# #     parameters$cd * parameters$drug_rel_c[drugs[successful]],
# #     treated_index
# #   )
# #   # Queue update to the last drug each successfully treated person received:
# #   variables$drug$queue_update(
# #     drugs[successful],
# #     treated_index
# #   )
# #   # Queue update to the time step on which each successfully treated person last
# #   # received a drug:
# #   variables$drug_time$queue_update(
# #     timestep,
# #     treated_index
# #   )
# # }
# 
# # Return the Bitset of treated individuals:
# treated_index
# 
# #----- 6a) Updated calculate_treated() function (V1): Resistance -> Efficacy -----------------------
# 
# calculate_treated_AM <- function(
#     variables,
#     clinical_infections,
#     parameters,
#     timestep,
#     renderer
# ) {
#   
#   # Gather the treatment coverage(s) in the current timestep:
#   treatment_coverages <- get_treatment_coverages(simparams, current_timestep)
#   
#   # Sum to get the total treatment coverage in the current time step:
#   ft <- sum(treatment_coverages)
#   
#   # If treatment = 0, return an empty Bitset for the treated individuals:
#   if (ft == 0) {
#     return(individual::Bitset$new(simparams$human_population))
#   }
#   
#   # Add the total treatment coverage to the renderer:
#   renderer$render('ft', ft, current_timestep)
#   
#   # Sample individuals from clinically infected to seek treatment using treatment coverage:
#   seek_treatment <- sample_bitset(clinical_infections, ft)
#   
#   # Store the number of people that seek treatment this time step:
#   n_treat <- seek_treatment$size()
#   
#   # Add the number of people seeking treatment in the current time step to the renderer:
#   renderer$render('n_treated', n_treat, current_timestep)
#   
#   # Assign each individual seeking treatment a drug based on their coverage(s):
#   drugs <- as.numeric(simparams$clinical_treatment_drugs[
#     sample.int(
#       length(simparams$clinical_treatment_drugs),
#       n_treat,
#       prob = treatment_coverages,
#       replace = TRUE
#     )
#   ])
#   
#   #+++ ANTIMALARIAL RESISTANCE +++#
#   
#   # Set the resistance drug index (each individuals drug):
#   AM_drug_index <- as.numeric(simparams$antimalarial_resistance_drug[drugs])
#   
#   # Open vectors to store the proportion of artemisinin resistance to the drug taken by each
#   # individual and it's probability of manifesting as early treatment failure:
#   art_resistance <- vector(); etf_prob <- vector()
#   
#   # For reach individual seeking treatment, retrieve the resistance and ETF probabilities corresponding
#   # to the drug they have been administered:
#   for(i in 1:n_treat) {
#     
#     # Identify the correct resistance index:
#     last_resistance <- match_timestep(ts = simparams$antimalarial_resistance_timesteps[[AM_drug_index[i]]],
#                                       t = current_timestep)
#     
#     # Retrieve the resistance value corresponding to the drug and time step:
#     art_resistance[i] <- simparams$prop_artemisinin_resistant[[AM_drug_index[i]]][last_resistance]
#     
#     # Retrieve the ETF phenotype probability corresponding to the drug and time step:
#     etf_prob[i] <- simparams$early_treatment_failure_prob[[AM_drug_index[i]]][last_resistance]
#     
#   }
#   
#   # Run bernoulli trials to determine which individuals experience early treatment failure given the
#   # drug they have been administered:
#   susceptible <- bernoulli_multi_p(p = 1 - (art_resistance * etf_prob))
#   
#   # Calculate the number of individuals who failed treatment due to early treatment failure:
#   n_ETF <- n_treat - length(susceptible)
#   
#   # Add the number of Early Treatment Failures to the renderer:
#   renderer$render('n_ETF', n_ETF, current_timestep)
#   
#   # Remove the individuals who failed treatment due to resistance from the drugs vector:
#   drugs <- drugs[susceptible] 
#   
#   # Subset the people who remained susceptible
#   susceptible_index <- bitset_at(seek_treatment, susceptible)
#   
#   #+++ DRUG EFFICACY +++#
#   # Determine, using drug efficacies, use bernoulli trials to see who's treatment worked:
#   successful <- bernoulli_multi_p(simparams$drug_efficacy[drugs])
#   
#   # Calculate the number of people who failed treatment due to efficacy:
#   n_treat_eff_fail <- length(susceptible) - length(successful)
#   
#   # Add the number of treatment efficacy failures to the renderer:
#   renderer$render('n_treat_eff_fail', n_treat_eff_fail, current_timestep)
#   
#   # Subset the successfully treated people:
#   treated_index <- bitset_at(susceptible_index, successful)
#   
#   # Add number of successfully treated individuals to the renderer
#   renderer$render('n_treat_success', treated_index$size(), current_timestep)
#   
#   # For those people who have been successfully treated:
#   if (treated_index$size() > 0) {
#     
#     # Queue update to the infectious state to Tr:
#     variables$state$queue_update('Tr', treated_index)
#     
#     # Queue update to the infectivity to reflect their new state (Tr):
#     variables$infectivity$queue_update(
#       parameters$cd * parameters$drug_rel_c[drugs[successful]],
#       treated_index
#     )
#     # Queue update to the last drug each successfully treated person received:
#     variables$drug$queue_update(
#       drugs[successful],
#       treated_index
#     )
#     # Queue update to the time step on which each successfully treated person last
#     # received a drug:
#     variables$drug_time$queue_update(
#       timestep,
#       treated_index
#     )
#   }
#   
#   # Return the Bitset of treated individuals:
#   treated_index
# }
# 
# #----- 6b) Updated calculate_treated() function (V1): Testing --------------------------------------
# 
# #----- X) Notes ------------------------------------------------------------------------------------
# 
# ##' If we are using a time-series of resistance proportions, we need a method for retrieving the up-
# ##' to-date proportions given the current time step and dates on which resistance values are updated.
# ##' For the time-varying deathrate, in create_mortality_process() (mortality_processes.R), the following
# ##' is used:
# ##'           last_deathrate <- match_timestep(parameters$deathrate_timesteps, timestep)
# ##'
# ##' We can use the same method in the following way:
# ##'
# ##'           latest_resistance <- match_timestep(unlist(parameters$antimalarial_resistance_timesteps),
# ##'                                               current_timestep)
# ##'
# 
# 
# 
# 
# 
# 

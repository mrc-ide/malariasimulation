test_that('set_antimalarial_resistance() toggles resistance on', {
  simparams <- get_parameters()
  simparams <- set_drugs(parameters = simparams, drugs = list(SP_AQ_params))
  simparams <- set_clinical_treatment(parameters = simparams,
                                      drug = 1,
                                      timesteps = 1,
                                      coverages = 1)
  set_antimalarial_resistance(parameters = simparams,
                              drug = 1,
                              timesteps = 1,
                              artemisinin_resistance = 0.5,
                              partner_drug_resistance = 0,
                              slow_parasite_clearance_prob = 0.5,
                              early_treatment_failure_prob = 0.6,
                              late_clinical_failure_prob = 0.2,
                              late_parasitological_prob = 0.3,
                              reinfection_prob = 0.4) -> simparams
  expect_identical(object = simparams$antimalarial_resistance, expected = TRUE)
})

test_that('set_antimalarial_resistance() errors if parameter inputs of different length to timesteps', {
  simparams <- get_parameters()
  simparams <- set_drugs(parameters = simparams, drugs = list(SP_AQ_params))
  simparams <- set_clinical_treatment(parameters = simparams,
                                      drug = 1,
                                      timesteps = 1,
                                      coverages = 1)
  expect_error(object =  set_antimalarial_resistance(parameters = simparams,
                                                     drug = 1,
                                                     timesteps = c(1, 10),
                                                     artemisinin_resistance = 0.5,
                                                     partner_drug_resistance = 0,
                                                     slow_parasite_clearance_prob = 0.5,
                                                     early_treatment_failure_prob = 0.6,
                                                     late_clinical_failure_prob = 0.2,
                                                     late_parasitological_prob = 0.3,
                                                     reinfection_prob = 0.4))
})

test_that('set_antimalarial_resistance() errors if resistance proportions outside bound of 0-1', {
  simparams <- get_parameters()
  simparams <- set_drugs(parameters = simparams, drugs = list(SP_AQ_params))
  simparams <- set_clinical_treatment(parameters = simparams,
                                      drug = 1,
                                      timesteps = 1,
                                      coverages = 1)
  expect_error(object =  set_antimalarial_resistance(parameters = simparams,
                                                     drug = 1,
                                                     timesteps = 1,
                                                     artemisinin_resistance = 1.01,
                                                     partner_drug_resistance = 0,
                                                     slow_parasite_clearance_prob = 0.5,
                                                     early_treatment_failure_prob = 0.6,
                                                     late_clinical_failure_prob = 0.2,
                                                     late_parasitological_prob = 0.3,
                                                     reinfection_prob = 0.4))
})

test_that('set_antimalarial_resistance() errors if resistance phenotype probabilities outside bound of 0-1', {
  simparams <- get_parameters()
  simparams <- set_drugs(parameters = simparams, drugs = list(SP_AQ_params))
  simparams <- set_clinical_treatment(parameters = simparams,
                                      drug = 1,
                                      timesteps = 1,
                                      coverages = 1)
  expect_error(object =  set_antimalarial_resistance(parameters = simparams,
                                                     drug = 1,
                                                     timesteps = 1,
                                                     artemisinin_resistance = 0.4,
                                                     partner_drug_resistance = 0,
                                                     slow_parasite_clearance_prob = -0.5,
                                                     early_treatment_failure_prob = 0.6,
                                                     late_clinical_failure_prob = 0.2,
                                                     late_parasitological_prob = 0.3,
                                                     reinfection_prob = 0.4))
})

test_that('set_antimalarial_resistance() errors if drug index > than number of drugs assigned using set_drugs()', {
  simparams <- get_parameters()
  simparams <- set_drugs(parameters = simparams, drugs = list(SP_AQ_params))
  simparams <- set_clinical_treatment(parameters = simparams,
                                      drug = 1,
                                      timesteps = 1,
                                      coverages = 1)
  expect_error(object =  set_antimalarial_resistance(parameters = simparams,
                                                     drug = 2,
                                                     timesteps = 1,
                                                     artemisinin_resistance = 0.4,
                                                     partner_drug_resistance = 0.3,
                                                     slow_parasite_clearance_prob = 0.5,
                                                     early_treatment_failure_prob = 0.6,
                                                     late_clinical_failure_prob = 0.2,
                                                     late_parasitological_prob = 0.3,
                                                     reinfection_prob = 0.4))
})

test_that('set_antimalarial_resistance() assigns parameters correctly despite drug order', {
  
  # Randomly assign the order of the drug assignment in clinical treatment and resistance:
  clinical_treatment_drugs <- sample.int(n = 3, size = 3, replace = FALSE)
  antimalarial_resistance_drugs <- sample.int(n = 3, size = 3, replace = FALSE)
  
  # Randomly assign drug coverages for clinical treatment:
  clinical_treatment_coverages <- sample(c(0.1, 0.2, 0.4))
  
  # Randomly assign timesteps for clinical treatment and resistance:
  clinical_treatment_timesteps <- sample.int(n = 100, size = 3, replace = FALSE)
  antimalarial_resistance_timesteps <- sample.int(n = 100, size = 3, replace = FALSE)
  
  # Randomly assign values to antimalarial resistance proportions
  antimalarial_resistance_artemisinin_proportion <- sample(x = seq(0, 0.33, by = 0.01), size = 3, replace = FALSE)
  antimalarial_resistance_partner_drug_proportion <- sample(x = seq(0, 0.33, by = 0.01), size = 3, replace = FALSE)
  
  # Randomly assign values to antimalarial resistance probabilities:
  antimalarial_resistance_ETF_probability <- sample(x = seq(0, 0.33, by = 0.01), size = 3, replace = FALSE)
  antimalarial_resistance_SPC_probability <- sample(x = seq(0, 0.33, by = 0.01), size = 3, replace = FALSE)
  antimalarial_resistance_LCF_probability <- sample(x = seq(0, 0.33, by = 0.01), size = 3, replace = FALSE)
  antimalarial_resistance_LPF_probability <- sample(x = seq(0, 0.33, by = 0.01), size = 3, replace = FALSE)
  antimalarial_resistance_RDP_probability <- sample(x = seq(0, 0.33, by = 0.01), size = 3, replace = FALSE)
  
  # Establish some parameters
  simparams <- get_parameters()
  
  # Add the three drug parameter vectors to the parameter list:
  simparams <- set_drugs(parameters = simparams, drugs = list(AL_params, SP_AQ_params, DHA_PQP_params))
  
  # Loop through and assign the parameters for each drug
  for(i in seq(clinical_treatment_drugs)) {
    simparams <- set_clinical_treatment(parameters = simparams,
                                        drug = clinical_treatment_drugs[i],
                                        timesteps = clinical_treatment_timesteps[i],
                                        coverages = clinical_treatment_coverages[i])
  }
  
  # Loop through and assign the antimalarial resistance parameters for each drug:
  for(i in seq(antimalarial_resistance_drugs)) {
    simparams <- set_antimalarial_resistance(parameters = simparams,
                                             drug = antimalarial_resistance_drugs[i],
                                             timesteps = antimalarial_resistance_timesteps[i],
                                             artemisinin_resistance = antimalarial_resistance_artemisinin_proportion[i],
                                             partner_drug_resistance = antimalarial_resistance_partner_drug_proportion[i],
                                             slow_parasite_clearance_prob = antimalarial_resistance_SPC_probability[i],
                                             early_treatment_failure_prob = antimalarial_resistance_ETF_probability[i],
                                             late_clinical_failure_prob = antimalarial_resistance_LCF_probability[i],
                                             late_parasitological_prob = antimalarial_resistance_LPF_probability[i],
                                             reinfection_prob = antimalarial_resistance_RDP_probability[i])
  }
  
  # Check the antimalarial resistance parameters after assignment using set_antimalarial_resistance():
  expect_identical(simparams$antimalarial_resistance, TRUE)
  expect_identical(unlist(simparams$antimalarial_resistance_drug), antimalarial_resistance_drugs)
  expect_identical(unlist(simparams$antimalarial_resistance_timesteps), antimalarial_resistance_timesteps)
  expect_identical(unlist(simparams$prop_artemisinin_resistant), antimalarial_resistance_artemisinin_proportion)
  expect_identical(unlist(simparams$prop_partner_drug_resistant), antimalarial_resistance_partner_drug_proportion)
  expect_identical(unlist(simparams$slow_parasite_clearance_prob), antimalarial_resistance_SPC_probability)
  expect_identical(unlist(simparams$early_treatment_failure_prob), antimalarial_resistance_ETF_probability)
  expect_identical(unlist(simparams$late_clinical_failure_prob), antimalarial_resistance_LCF_probability)
  expect_identical(unlist(simparams$late_parasitological_failure_prob), antimalarial_resistance_LPF_probability)
  expect_identical(unlist(simparams$reinfection_during_prophylaxis), antimalarial_resistance_RDP_probability)
  
})

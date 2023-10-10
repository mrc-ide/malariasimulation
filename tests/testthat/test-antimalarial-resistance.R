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

test_that('set_antimalarial_resistance() errors if resistance proportions outside of range 0-1', {
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
                                                     reinfection_prob = 0.4), 
               regexp = "Artemisinin and partner-drug resistance proportions must fall between 0 and 1")
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

test_that('set_antimalarial_resistance() assigns parameters correctly despite order in which resistance parameters are specified', {
  
  parameters <- get_parameters()
  parameters <- set_drugs(parameters = parameters, drugs = list(AL_params, SP_AQ_params, DHA_PQP_params))
  parameters <- set_clinical_treatment(parameters = parameters, drug = 2, timesteps = 1, coverages = 0.2)
  parameters <- set_clinical_treatment(parameters = parameters, drug = 1, timesteps = 1, coverages = 0.1)
  parameters <- set_clinical_treatment(parameters = parameters, drug = 3, timesteps = 1, coverages = 0.4)
  parameters <- set_antimalarial_resistance(parameters = parameters,
                                            drug = 2,
                                            timesteps = 1,
                                            artemisinin_resistance = 0.5,
                                            partner_drug_resistance = 0,
                                            slow_parasite_clearance_prob = 0.41,
                                            early_treatment_failure_prob = 0.2,
                                            late_clinical_failure_prob = 0,
                                            late_parasitological_prob = 0,
                                            reinfection_prob = 0)
  parameters <- set_antimalarial_resistance(parameters = parameters,
                                            drug = 3,
                                            timesteps = 1,
                                            artemisinin_resistance = 0,
                                            partner_drug_resistance = 0.43,
                                            slow_parasite_clearance_prob = 0,
                                            early_treatment_failure_prob = 0,
                                            late_clinical_failure_prob = 0.01,
                                            late_parasitological_prob = 0.42,
                                            reinfection_prob = 0.89)
  parameters <- set_antimalarial_resistance(parameters = parameters,
                                            drug = 1,
                                            timesteps = 1,
                                            artemisinin_resistance = 0.27,
                                            partner_drug_resistance = 0.61,
                                            slow_parasite_clearance_prob = 0.23,
                                            early_treatment_failure_prob = 0.9,
                                            late_clinical_failure_prob = 0.49,
                                            late_parasitological_prob = 0.81,
                                            reinfection_prob = 0.009)
  
  expect_identical(parameters$antimalarial_resistance, TRUE)
  expect_identical(unlist(parameters$antimalarial_resistance_drug), c(2, 3, 1))
  expect_identical(unlist(parameters$antimalarial_resistance_timesteps), rep(1, 3))
  expect_identical(unlist(parameters$prop_artemisinin_resistant), c(0.5, 0, 0.27))
  expect_identical(unlist(parameters$prop_partner_drug_resistant), c(0, 0.43, 0.61))
  expect_identical(unlist(parameters$slow_parasite_clearance_prob), c(0.41, 0, 0.23))
  expect_identical(unlist(parameters$early_treatment_failure_prob), c(0.2, 0, 0.9))
  expect_identical(unlist(parameters$late_clinical_failure_prob), c(0, 0.01, 0.49))
  expect_identical(unlist(parameters$late_parasitological_failure_prob), c(0, 0.42, 0.81))
  expect_identical(unlist(parameters$reinfection_during_prophylaxis), c(0, 0.89, 0.009))
  
})

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
                              late_clinical_failure_prob = 0,
                              late_parasitological_prob = 0,
                              reinfection_prob = 0, 
                              slow_parasite_clearance_time = 10) -> simparams
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
                                                     reinfection_prob = 0.4,
                                                     slow_parasite_clearance_time = 10))
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
                                                     late_clinical_failure_prob = 0,
                                                     late_parasitological_prob = 0,
                                                     reinfection_prob = 0,
                                                     slow_parasite_clearance_time = 10), 
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
                                                     reinfection_prob = 0.4, 
                                                     slow_parasite_clearance_time = 5))
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
                                            reinfection_prob = 0, 
                                            slow_parasite_clearance_time = 5)
  parameters <- set_antimalarial_resistance(parameters = parameters,
                                            drug = 3,
                                            timesteps = 1,
                                            artemisinin_resistance = 0,
                                            partner_drug_resistance = 0,
                                            slow_parasite_clearance_prob = 0,
                                            early_treatment_failure_prob = 0,
                                            late_clinical_failure_prob = 0,
                                            late_parasitological_prob = 0,
                                            reinfection_prob = 0, 
                                            slow_parasite_clearance_time = 10)
  parameters <- set_antimalarial_resistance(parameters = parameters,
                                            drug = 1,
                                            timesteps = 1,
                                            artemisinin_resistance = 0.27,
                                            partner_drug_resistance = 0,
                                            slow_parasite_clearance_prob = 0.23,
                                            early_treatment_failure_prob = 0.9,
                                            late_clinical_failure_prob = 0,
                                            late_parasitological_prob = 0,
                                            reinfection_prob = 0, 
                                            slow_parasite_clearance_time = 20)
  
  expect_identical(parameters$antimalarial_resistance, TRUE)
  expect_identical(unlist(parameters$antimalarial_resistance_drug), c(2, 3, 1))
  expect_identical(unlist(parameters$antimalarial_resistance_timesteps), rep(1, 3))
  expect_identical(unlist(parameters$prop_artemisinin_resistant), c(0.5, 0, 0.27))
  expect_identical(unlist(parameters$prop_partner_drug_resistant), c(0, 0, 0))
  expect_identical(unlist(parameters$slow_parasite_clearance_prob), c(0.41, 0, 0.23))
  expect_identical(unlist(parameters$early_treatment_failure_prob), c(0.2, 0, 0.9))
  expect_identical(unlist(parameters$late_clinical_failure_prob), c(0, 0, 0))
  expect_identical(unlist(parameters$late_parasitological_failure_prob), c(0, 0, 0))
  expect_identical(unlist(parameters$reinfection_during_prophylaxis), c(0, 0, 0))
  expect_identical(unlist(parameters$dt_slow_parasite_clearance), c(5, 10, 20))
  
})

test_that(desc = "set_antimalarial_resistance errors if length slow_parasite_clearance_time > 1", code = {
  
  parameters <- get_parameters()
  parameters <- set_drugs(parameters = parameters, drugs = list(SP_AQ_params))
  parameters <- set_clinical_treatment(parameters = parameters,
                                       drug = 1,
                                       timesteps = c(0, 10), 
                                       coverages = c(0.1, 0.2))
  
  expect_error(
    parameters <- set_antimalarial_resistance(parameters = parameters,
                                              drug = 1,
                                              timesteps = c(0, 10),
                                              artemisinin_resistance = c(0.4, 0.8),
                                              partner_drug_resistance = c(0, 0),
                                              slow_parasite_clearance_prob = c(0.2, 0.4),
                                              early_treatment_failure_prob = c(0, 0.45),
                                              late_clinical_failure_prob = c(0, 0),
                                              late_parasitological_prob = c(0, 0),
                                              reinfection_prob = c(0, 0), 
                                              slow_parasite_clearance_time = c(10 ,11)),
    "Error: length of slow_parasite_clearance_time not equal to 1")
})

test_that(desc = "set_antimalarial_resistance errors if slow_parasite_clearance_time not positive", code = {
  
  parameters <- get_parameters()
  parameters <- set_drugs(parameters = parameters, drugs = list(SP_AQ_params))
  parameters <- set_clinical_treatment(parameters = parameters,
                                       drug = 1,
                                       timesteps = c(0, 10), 
                                       coverages = c(0.1, 0.2))
  
  expect_error(
    parameters <- set_antimalarial_resistance(parameters = parameters,
                                              drug = 1,
                                              timesteps = c(0, 10),
                                              artemisinin_resistance = c(0.4, 0.8),
                                              partner_drug_resistance = c(0, 0),
                                              slow_parasite_clearance_prob = c(0.2, 0.4),
                                              early_treatment_failure_prob = c(0, 0.45),
                                              late_clinical_failure_prob = c(0, 0),
                                              late_parasitological_prob = c(0, 0),
                                              reinfection_prob = c(0, 0), 
                                              slow_parasite_clearance_time = c(0)),
    "Error: slow_parasite_clearance_time is non-positive")
})

test_that('get_antimalarial_resistance_parameters() correctly retrieves parameters when multiple drugs assigned', {
    
    get_parameters(overrides = list(human_population = 10000)) %>%
      set_drugs(drugs = list(AL_params, SP_AQ_params, DHA_PQP_params)) %>%
      set_clinical_treatment(drug = 1, timesteps = 1, coverages = 0.4) %>%
      set_clinical_treatment(drug = 2, timesteps = 1, coverages = 0.3) %>%
      set_clinical_treatment(drug = 3, timesteps = 1, coverages = 0.2) %>%
      set_equilibrium(init_EIR = 20) %>%
      set_antimalarial_resistance(drug = 2,
                                  timesteps = c(0, 20),
                                  artemisinin_resistance = c(0.02, 0.2),
                                  partner_drug_resistance = c(0, 0),
                                  slow_parasite_clearance_prob = c(0.02, 0.2),
                                  early_treatment_failure_prob = c(0.02, 0.2),
                                  late_clinical_failure_prob = c(0, 0),
                                  late_parasitological_prob = c(0, 0),
                                  reinfection_prob = c(0, 0), 
                                  slow_parasite_clearance_time = 20) %>%
      set_antimalarial_resistance(drug = 1,
                                  timesteps = c(0, 10),
                                  artemisinin_resistance = c(0.01, 0.1),
                                  partner_drug_resistance = c(0, 0),
                                  slow_parasite_clearance_prob = c(0.01, 0.1),
                                  early_treatment_failure_prob = c(0.01, 0.1),
                                  late_clinical_failure_prob = c(0, 0),
                                  late_parasitological_prob = c(0, 0),
                                  reinfection_prob = c(0, 0),
                                  slow_parasite_clearance_time = 10) %>%
      set_antimalarial_resistance(drug = 3,
                                  timesteps = c(0, 30),
                                  artemisinin_resistance = c(0.03, 0.3),
                                  partner_drug_resistance = c(0, 0),
                                  slow_parasite_clearance_prob = c(0.03, 0.3),
                                  early_treatment_failure_prob = c(0.03, 0.3),
                                  late_clinical_failure_prob = c(0, 0),
                                  late_parasitological_prob = c(0, 0),
                                  reinfection_prob = c(0, 0),
                                  slow_parasite_clearance_time = 30) -> parameters
    
    drugs <- c(1, 3, 2, 1, 2, 3, 3, 3, 2, 1, 3, 1, 2, 3, 2)
    timestep <- 25
    
    resistance_parameters <- get_antimalarial_resistance_parameters(parameters = parameters,
                                                                    drugs = drugs, 
                                                                    timestep = timestep)
    
    expected_resistance_parameters <- list()
    expected_resistance_parameters$artemisinin_resistance_proportion <- c(0.1, 0.03, 0.2, 0.1, 0.2, 0.03, 0.03, 0.03, 0.2, 0.1, 0.03, 0.1, 0.2, 0.03, 0.2)
    expected_resistance_parameters$partner_drug_resistance_proportion <- rep(0, 15)
    expected_resistance_parameters$slow_parasite_clearance_probability <- c(0.1, 0.03, 0.2, 0.1, 0.2, 0.03, 0.03, 0.03, 0.2, 0.1, 0.03, 0.1, 0.2, 0.03, 0.2)
    expected_resistance_parameters$early_treatment_failure_probability <- c(0.1, 0.03, 0.2, 0.1, 0.2, 0.03, 0.03, 0.03, 0.2, 0.1, 0.03, 0.1, 0.2, 0.03, 0.2)
    expected_resistance_parameters$late_clinical_failure_probability <- rep(0, 15)
    expected_resistance_parameters$late_parasitological_failure_probability <- rep(0, 15)
    expected_resistance_parameters$reinfection_during_prophylaxis_probability <- rep(0, 15)
    expected_resistance_parameters$dt_slow_parasite_clearance <- c(10, 30, 20, 10, 20, 30, 30, 30, 20, 10, 30, 10, 20, 30, 20)
    
    expect_identical(resistance_parameters, expected = expected_resistance_parameters)
    
  })

test_that('get_antimalarial_resistance_parameters() correctly retrieves parameters when not all drugs assigned resistance', {
    
    get_parameters(overrides = list(human_population = 10000)) %>%
      set_drugs(drugs = list(AL_params, SP_AQ_params, DHA_PQP_params)) %>%
      set_clinical_treatment(drug = 1, timesteps = 1, coverages = 0.4) %>%
      set_clinical_treatment(drug = 2, timesteps = 1, coverages = 0.3) %>%
      set_clinical_treatment(drug = 3, timesteps = 1, coverages = 0.2) %>%
      set_equilibrium(init_EIR = 20) %>%
      set_antimalarial_resistance(drug = 2,
                                  timesteps = c(0, 20),
                                  artemisinin_resistance = c(0.02, 0.2),
                                  partner_drug_resistance = c(0, 0),
                                  slow_parasite_clearance_prob = c(0.02, 0.2),
                                  early_treatment_failure_prob = c(0.02, 0.2),
                                  late_clinical_failure_prob = c(0, 0),
                                  late_parasitological_prob = c(0, 0),
                                  reinfection_prob = c(0, 0), 
                                  slow_parasite_clearance_time = 20) -> parameters
    
    drugs <- c(1, 3, 2, 1, 2, 3, 3, 3, 2, 1, 3, 1, 2, 3, 2)
    timestep <- 25

        resistance_parameters <- get_antimalarial_resistance_parameters(parameters = parameters,
                                                                    drugs = drugs, 
                                                                    timestep = timestep)
    
    expected_resistance_parameters <- list()
    expected_resistance_parameters$artemisinin_resistance_proportion <- c(0, 0, 0.2, 0, 0.2, 0, 0, 0, 0.2, 0, 0, 0, 0.2, 0, 0.2)
    expected_resistance_parameters$partner_drug_resistance_proportion <- rep(0, 15)
    expected_resistance_parameters$slow_parasite_clearance_probability <- c(0, 0, 0.2, 0, 0.2, 0, 0, 0, 0.2, 0, 0, 0, 0.2, 0, 0.2)
    expected_resistance_parameters$early_treatment_failure_probability <- c(0, 0, 0.2, 0, 0.2, 0, 0, 0, 0.2, 0, 0, 0, 0.2, 0, 0.2)
    expected_resistance_parameters$late_clinical_failure_probability <- rep(0, 15)
    expected_resistance_parameters$late_parasitological_failure_probability <- rep(0, 15)
    expected_resistance_parameters$reinfection_during_prophylaxis_probability <- rep(0, 15)
    expected_resistance_parameters$dt_slow_parasite_clearance <- c(5, 5, 20, 5, 20, 5, 5, 5, 20, 5, 5, 5, 20, 5, 20)
    
    expect_identical(resistance_parameters, expected = expected_resistance_parameters)
    
  })

test_that('get_antimalarial_resistance_parameters() returns an error when antimalarial resistance has not been parameterised', {
    
    get_parameters(overrides = list(human_population = 10000)) %>%
      set_drugs(drugs = list(AL_params, SP_AQ_params, DHA_PQP_params)) %>%
      set_clinical_treatment(drug = 1, timesteps = 1, coverages = 0.4) %>%
      set_clinical_treatment(drug = 2, timesteps = 1, coverages = 0.3) %>%
      set_clinical_treatment(drug = 3, timesteps = 1, coverages = 0.2) %>%
      set_equilibrium(init_EIR = 20) -> parameters
    
    drugs <- c(1, 3, 2, 1, 2, 3, 3, 3, 2, 1, 3, 1, 2, 3, 2)
    timestep <- 25
    
    
    expect_error(get_antimalarial_resistance_parameters(parameters = parameters,
                                                        drugs = drugs,
                                                        timestep = timestep),
                 "Error: Antimalarial resistance has not been parameterised; antimalarial_resistance = FALSE")
})
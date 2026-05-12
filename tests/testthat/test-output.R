test_that('Age-stratified state counts are present in output', {
  parameters <- get_parameters()
  parameters <- set_epi_outputs(parameters, state_count = c(0, 1825))
  sim <- run_simulation(timesteps = 5, parameters)

  states <- c('S', 'A', 'D', 'U', 'Tr')
  expected_cols <- paste0(states, '_count_0_1824')
  expect_in(expected_cols, names(sim))
})

test_that('Age-stratified state counts sum to n_age for that group', {
  parameters <- get_parameters()
  parameters <- set_epi_outputs(
    parameters,
    age_group = c(0, 1825),
    state_count = c(0, 1825)
  )
  sim <- run_simulation(timesteps = 5, parameters)

  state_total <- sim$S_count_0_1824 + sim$A_count_0_1824 + sim$D_count_0_1824 +
    sim$U_count_0_1824 + sim$Tr_count_0_1824
  expect_equal(state_total, sim$n_age_0_1824)
})

test_that('Age-stratified PCR-detectable count matches state counts', {
  parameters <- get_parameters()
  parameters <- set_epi_outputs(
    parameters,
    prevalence = c(0, 1825),
    state_count = c(0, 1825)
  )
  sim <- run_simulation(timesteps = 5, parameters)

  pcr_from_states <- sim$A_count_0_1824 + sim$D_count_0_1824 +
    sim$U_count_0_1824 + sim$Tr_count_0_1824
  expect_equal(pcr_from_states, sim$n_detect_pcr_0_1824)
})

test_that('Age-stratified state counts respect multiple age groups', {
  parameters <- get_parameters()
  parameters <- set_epi_outputs(
    parameters,
    age_group = c(0, 1825, 3650),
    state_count = c(0, 1825, 3650)
  )
  sim <- run_simulation(timesteps = 5, parameters)

  # Each age group sums to its n_age
  for (bounds in list(c(0, 1824), c(1825, 3649))) {
    lower <- bounds[1]
    upper <- bounds[2]
    state_total <- sim[[paste0('S_count_', lower, '_', upper)]] +
      sim[[paste0('A_count_', lower, '_', upper)]] +
      sim[[paste0('D_count_', lower, '_', upper)]] +
      sim[[paste0('U_count_', lower, '_', upper)]] +
      sim[[paste0('Tr_count_', lower, '_', upper)]]
    expect_equal(state_total, sim[[paste0('n_age_', lower, '_', upper)]])
  }

  # Age groups sum to total population
  total_from_groups <- (sim$S_count_0_1824 + sim$A_count_0_1824 + sim$D_count_0_1824 +
    sim$U_count_0_1824 + sim$Tr_count_0_1824) +
    (sim$S_count_1825_3649 + sim$A_count_1825_3649 + sim$D_count_1825_3649 +
    sim$U_count_1825_3649 + sim$Tr_count_1825_3649)
  expect_true(all(total_from_groups <= sim$S_count + sim$A_count + sim$D_count +
    sim$U_count + sim$Tr_count))
})

test_that('Test age parameter function works', {
  parameters <- get_parameters()
  age_limits <- c(0,1,2,3)*365
  parameters <- set_epi_outputs(parameters,
                                age_group = age_limits,
                                incidence = age_limits+1,
                                clinical_incidence = age_limits+2,
                                severe_incidence = age_limits+3,
                                prevalence = age_limits+4,
                                ica = age_limits+5,
                                icm = age_limits+6,
                                id = age_limits+7,
                                ib = age_limits+8,
                                iva = age_limits+9,
                                ivm = age_limits+10
  )
  
  sim <- run_simulation(timesteps = 1, parameters)
  
  prefixes <- c("n_age",
                "n_inc", "p_inc",
                "n_inc_clinical","p_inc_clinical",
                "n_inc_severe","p_inc_severe",
                "n_detect_lm","p_detect_lm","n_detect_pcr",
                "ica_mean", "icm_mean","id_mean","ib_mean","iva_mean","ivm_mean")
  offsets <- c(0, rep(1, 2), rep(2, 2), rep(3, 2), rep(4, 3), 5:10)
  expect_equal(length(prefixes), length(offsets))

  expected <- paste0(rep(prefixes, each = 3),
                     "_", age_limits[-4]+rep(offsets, each = 3),
                     "_", age_limits[-1]-1+rep(offsets, each = 3))
  expect_in(expected, names(sim))
  
  
  expect_error(set_epi_outputs(parameters, iaa = age_limits))
  expect_error(get_parameters(parasite = "vivax") |> 
                 set_epi_outputs(iva = age_limits))
  
})

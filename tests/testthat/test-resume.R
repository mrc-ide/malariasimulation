test_that('Simulation can be resumed', {
  initial_timesteps <- 50
  total_timesteps <- 100

  parameters <- get_parameters()

  set.seed(1)
  first_phase <- run_resumable_simulation(initial_timesteps, parameters)
  second_phase <- run_resumable_simulation(
    total_timesteps,
    parameters,
    initial_state=first_phase$state,
    restore_random_state=TRUE)

  set.seed(1)
  uninterrupted_run <- run_simulation(total_timesteps, parameters=parameters)

  expect_equal(nrow(first_phase$data), initial_timesteps)
  expect_equal(nrow(second_phase$data), total_timesteps - initial_timesteps)
  expect_equal(rbind(first_phase$data, second_phase$data), uninterrupted_run)
})

test_that('Intervention parameters can be changed when resuming', {
  initial_timesteps <- 50
  total_timesteps <- 100
  tbv_timesteps <- 70

  # Because of how event scheduling works, we must enable TBV even in the inital phase.
  # We set a coverage of 0 to act as-if it was disabled.
  initial_parameters <- get_parameters() |> set_tbv(timesteps=tbv_timesteps, coverage=0, ages=5:60)

  tbv_parameters <- initial_parameters |>
    set_tbv(timesteps=tbv_timesteps, coverage=1, ages=5:60)

  initial_run <- run_resumable_simulation(initial_timesteps, initial_parameters)
  control_run <- run_resumable_simulation(total_timesteps, initial_parameters, initial_state = initial_run$state)
  tbv_run <- run_resumable_simulation(total_timesteps, tbv_parameters, initial_state = initial_run$state)

  expect_equal(control_run$data$n_vaccinated_tbv[tbv_timesteps - initial_timesteps], 0)
  expect_gt(tbv_run$data$n_vaccinated_tbv[tbv_timesteps - initial_timesteps], 0)
})

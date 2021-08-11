test_that('Adult ODE jumping doesnt get weirdly stuck', {
  parameters <- get_parameters(list(
    human_population = 10000,
    individual_mosquitoes = FALSE,
    model_seasonality = TRUE,
    g0 = 1.639449,
    g = c(-2.1978530, 0.6981486, 0.1225242),
    h = c(-1.4126293, 1.5910602, -0.7558194),
    average_age = 12480.69,
    Q0 = 0.01797597
  ))
  parameters <- set_equilibrium(parameters, 25.36424)
  models <- parameterise_mosquito_models(parameters)
  solver <- parameterise_solvers(models, parameters)[[1]]
  timestep <- 4314
  solver_jump(solver, timestep)
  solver_step(solver)
  expect_true(TRUE)
})

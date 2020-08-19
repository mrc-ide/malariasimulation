test_that('Seasonality correctly affects P', {
  print('helllloooo')
  parameters <- get_parameters(list(
    model_seasonality = TRUE,
    g0 = 0.285505,
    g = c(-0.325352, -0.0109352, 0.0779865),
    h = c(-0.132815, 0.104675, -0.013919),
    variety_proportions = 1
  ))
  total_M <- 1000
  model <- parameterise_ode(parameters)[[1]]
  timesteps <- 365 * 1

  counts <- c()
  
  for (t in seq(timesteps)) {
    print(paste0('stepping ', t))
    counts <- rbind(counts, c(t, mosquito_model_get_states(model)))
    mosquito_model_step(model, total_M)
  }

  df <- as.data.frame(counts)
  names(df) <- c('timestep', 'E', 'L', 'P')
  df$K <- vnapply(
    df$timestep,
    function(t) {
      carrying_capacity(
        t,
        TRUE,
        1,
        parameters$g0,
        parameters$g,
        parameters$h,
        calculate_carrying_capacity(parameters),
        calculate_R_bar(parameters)
      )
    }
  )
  ggplot2::ggplot(df) + ggplot2::geom_line(ggplot2::aes(x = timestep, y = P))
  ggplot2::ggsave('Pv.png')

  ggplot2::ggplot(df) + ggplot2::geom_line(ggplot2::aes(x = timestep, y = K))
  ggplot2::ggsave('K.png')
})

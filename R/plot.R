plot_counts_against_equilibrium <- function(output, pop_size, state_props) {
  sim_length <- dim(output)[[1]]
  equilibrium <- data.frame(
    timestep = seq(sim_length),
    human_S_count = rep(state_props[['S']] * pop_size, sim_length),
    human_D_count = rep(state_props[['D']] * pop_size, sim_length),
    human_A_count = rep(state_props[['A']] * pop_size, sim_length),
    human_U_count = rep(state_props[['U']] * pop_size, sim_length)
  )
  if ('repetitions' %in% colnames(output)) {
    repetitions <- max(output$repetitions)
  } else {
    repetitions <- 1
  }
  ggplot(
    melt(output[c(
      'timestep',
      'human_S_count',
      'human_D_count',
      'human_A_count',
      'human_U_count'
    )],
    'timestep'
    ),
  ) + geom_line(
    aes(x = timestep, y = value, group = variable, color = variable, alpha = 1/repetitions)
  ) + geom_line(
    data = melt(equilibrium, 'timestep'),
    aes(x = timestep, y = value, group = variable, color = variable),
    linetype = 'dashed'
  )
}

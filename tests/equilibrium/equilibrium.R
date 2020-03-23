library(malariasimulation)
library(malariaEquilibrium)
library(ggplot2)

remove_keys <- function(x, n) { for (name in n) { x[[name]] <- NULL }; x }

jamie_params <- load_parameter_set("Jamie_parameters.rds")

params <- remove_keys(
  jamie_params,
  c(
    's2',
    'rT',
    'rP',
    'tl',
    'g_inf',
    'fd0',
    'aA',
    'aU',
    'b1',
    'PM',
    'tau',
    'mu',
    'f',
    'Q0',
    'cd_w',
    'cd_p',
    'cT'
  )
)

simparams <- translate_jamie(params)

output <- run_simulation(300, simparams)

plot_states <- function(output) {
    # group the state counts into one column
    cols <- colnames(output)[colnames(output) != 'timestep']
    reshaped <- reshape(
        output,
        varying = cols,
        idvar = 'timestep',
        timevar = 'state',
        direction = 'long',
        v.names = c('counts')
    )

    # make the state column more readable
    reshaped$state <- vapply(reshaped$state, function(i) cols[[i]], character(1))

    ggplot(
        reshaped,
        aes(x = timestep, y = counts, group = state)
    ) + geom_line(aes(color = state))
}

plot_states(output[c(
  'human_A_count',
  'human_D_count',
  'human_S_count',
  'human_U_count',
  'human_I_count'
)])

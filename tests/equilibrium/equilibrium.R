library(malariasimulation)
library(malariaEquilibrium)
library(ggplot2)

sim_length <- 1000
burn_in <- 500

remove_keys <- function(x, n) { for (name in n) { x[[name]] <- NULL }; x }

jamie_params <- load_parameter_set("Jamie_parameters.rds")

params <- remove_keys(
  jamie_params,
  c(
    's2',
    'rT', # makes sense
    'rP', # makes sense
    'tl',
    'g_inf',
    'fd0',
    'aA',
    'aU',
    'b1',
    'PM',
    'tau',
    #'mu', mum
    'f',
    'Q0',
    'cd_w',
    'cd_p',
    'cT',
    'dE' # not sure if this translation works
  )
)

simparams <- translate_jamie(params)

output <- run_simulation(sim_length, simparams)

# Estimating EIR of the model
# Leave a 100 timestep grace period for the EIR to flatten out
ggplot(
  output,
  aes(x = timestep, y = EIR)
) + geom_line()

EIR <- mean(output$EIR[output$timestep > burn_in])

print(paste("Estimated equilibrium to be", EIR, sep=" "))

# Calculate equilibrium
eq <- human_equilibrium(EIR = EIR, ft = 0, p = jamie_params, age = 0:100)
state_props <- colSums(eq$states[,c('S', 'D', 'A', 'U')])

simparams$s_proportion <- state_props[['S']]
simparams$d_proportion <- state_props[['D']]
simparams$a_proportion <- state_props[['A']]
simparams$u_proportion <- state_props[['U']]

# Run from equilibrium
output <- run_simulation(sim_length, simparams)

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
  'human_U_count'
)])

plot_states(output[c(
  'mosquito_E_count',
  'mosquito_L_count',
  'mosquito_P_count',
  'mosquito_Im_count',
  'mosquito_Sm_count'
)])
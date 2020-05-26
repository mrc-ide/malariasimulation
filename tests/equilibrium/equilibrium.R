library(malariasimulation)
library(malariaEquilibrium)
library(ggplot2)
library(reshape2)

sim_length <- 10 * 365

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
    'f',
    'Q0',
    'cd_w',
    'cd_p',
    'cT',
    'dE' # not sure if this translation works
  )
)

simparams <- translate_jamie(params)

#add blood meal rates
simparams[c('av1', 'av2', 'av3')] <- jamie_params$f * jamie_params$Q0

output <- run_simulation(sim_length, simparams)

# Estimating EIR of the model
# Leave a burn in for the EIR to flatten out
ggplot(
  subset(output, output$timestep > (sim_length - 365)),
  aes(x = timestep, y = mean_EIR)
) + geom_line()

set.seed(42)
EIR <- mean(output$mean_EIR[output$timestep > (sim_length - 365)])

print(paste("Estimated equilibrium to be", EIR, sep=" "))

# Calculate equilibrium
eq <- human_equilibrium(EIR = EIR*365, ft = 0, p = jamie_params, age = 0:100)
state_props <- colSums(eq$states[,c('S', 'D', 'A', 'U')])

simparams$s_proportion <- state_props[['S']]
simparams$d_proportion <- state_props[['D']]
simparams$a_proportion <- state_props[['A']]
simparams$u_proportion <- state_props[['U']]

# Run from equilibrium
set.seed(42)
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

plot_states(subset(output, output$timestep > (sim_length - 365))[c(
  'mosquito_Im_count',
  'mosquito_Sm_count'
)])

ggplot(
  melt(
    output[c(
      'timestep',
      'human_ICA_mean',
      'human_ICM_mean',
      'human_IB_mean'
      )],
    'timestep'
  ),
  aes(x = timestep, y = value, color = variable)
) + geom_line()

EIR <- mean(output$mean_EIR[output$timestep > (sim_length - 365)])
print(paste("Actual EIR:", EIR, sep=" "))

library(malariasimulation)
library(malariaEquilibrium)

sim_length <- 1000

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

output <- run_simulation(sim_length, simparams)

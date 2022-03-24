
# use malariaEquilibrium
# (documented here: https://github.com/mrc-ide/malariaEquilibrium)
eirs <- seq(from = 10, to = 40, by=.1)
eq_params <- malariaEquilibrium::load_parameter_set("Jamie_parameters.rds")

# calculate prevalence between 2:10 for a bunch of EIRs
prevs <- vapply(
  eirs,
  function(eir) {
    eq <- malariaEquilibrium::human_equilibrium(
      eir,
      ft=0,
      p=eq_params,
      age=0:100
    )
    sum(eq$states[2:10, 'pos_M']) / sum(eq$states[2:10, 'prop'])
  },
  numeric(1)
)

# have a look at the eir prev relationship
plot(eirs, prevs)

# calculate total_M for a bunch of EIRs
total_Ms <- vapply(
  eirs,
  function(eir) {
    malariasimulation::set_equilibrium(
      malariasimulation::get_parameters(),
      eir,
      eq_params=eq_params
    )$total_M
  },
  numeric(1)
)

# have a look at the eir total_M relationship
plot(eirs, total_Ms)


# set up your simulation with the EIR you want
year <- 365
EIR <- 16

# get default parameters
params <- malariasimulation::get_parameters()
params <- malariasimulation::set_equilibrium(params, EIR, eq_params=eq_params)

output <- malariasimulation::run_simulation(year * 5, parameters = params)

# plot some outputs
ggplot(output) + geom_line(aes(timestep, pv_730_3650))
ggplot(output) + geom_line(aes(timestep, total_M))
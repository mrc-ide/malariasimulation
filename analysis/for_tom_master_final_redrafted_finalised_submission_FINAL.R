#--------------------------------------------------------------------------------------------------#
#----- for_tom_master_final_redrafted_finalised_submission.R -----#
#-----------------------------------------------------------------#

# Load the requisite packages:
library(malariasimulation)
library(tidyverse)
library(cowplot)

#----- 1) Simulation Set-Up ------------------------------------------------------------------------

# Set up basic simulation parameters list:
year <- 365
sim_length <- 365*10
human_population <- 1000
starting_EIR <- 150
simparams <- get_parameters(overrides = list(
  human_population = human_population,
  prevalence_rendering_min_ages = 0,
  prevalence_rendering_max_ages = 5 * 365,
  individual_mosquitoes = FALSE
))

# Mosquito parameters:
mosq_params <- gamb_params
mosq_params$Q0 <- 0.9
mosq_params$phi_bednets <- 0.89

# Update parameter list with mosquito parameters:
simparams <- set_species(
  simparams,
  species = list(mosq_params),
  proportions = c(1)
)

# Calibrate the initial conditions to the parameter set:
simparams <- set_equilibrium(
  parameters = simparams,
  init_EIR = starting_EIR
)

# Parameterise the simulations for a bed-net campaign:
bednetstimesteps <- c(1, 4, 7, 10)*365
bednetparams <- set_bednets(
  simparams,
  timesteps = bednetstimesteps,
  coverages = c(0.75, 0.75, 0.75, 0.75),
  retention = 1000*365,
  dn0 = matrix(c(0.41, 0.41, 0.41, 0.41), nrow = 4, ncol = 1),
  rn = matrix(c(0.56, 0.56, 0.56, 0.56), nrow = 4, ncol = 1),
  rnm = matrix(c(0.24, 0.24, 0.24, 0.24), nrow = 4, ncol = 1),
  gamman = rep(2.64*365, 4)
)

# Set up correlation parameters to govern who receives the bednets:
correlations <- get_correlation_parameters(
  parameters = bednetparams
)
correlations$inter_round_rho('bednets', 1)

#----- 2) Version with rewritten run_simulations_with_repetitions() function -----------------------

##' The updated version of run_simulation_with_repetitions() replaces the overrides argument with
##' parameters and correlations arguments, which are now fed directly to the internal run_simulations()
##' function call.

# Store the updated version of run_simulation_with_repetitions():
updated_run_simulation_with_repetitions <- function(
    timesteps,
    repetitions,
    parameters,
    correlations,
    parallel = FALSE
) {
  if (parallel) {
    fapply <- parallel::mclapply
  } else {
    fapply <- lapply
  }
  dfs <- fapply(
    seq(repetitions),
    function(repetition) {
      df <- run_simulation(
        timesteps = timesteps,
        parameters = parameters,
        correlations = correlations
      )
      df$repetition <- repetition
      df
    }
  )
  do.call("rbind", dfs)
}

# Run the simulations using the updated_run_simulation_with_repetitions() function:
outputs_v1 <- updated_run_simulation_with_repetitions(
  timesteps = sim_length,
  repetitions = 5,
  parameters = bednetparams,
  correlations = correlations,
  parallel = FALSE
)

# Plot the PfPR[2-10] through time for each iteration:
outputs_v1 |>
  mutate(pfpr = n_detect_lm_0_1825 / n_age_0_1825) |>
  ggplot(aes(x = timestep, y = pfpr, colour = as.factor(repetition))) + geom_line() +
  theme_bw() +
  labs(x = "Time (days)", y = "Prevalence", colour = "Iteration") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

#----- 3) Version with custom run_simulation loop function -----------------------------------------

timesteps <- sim_length
iterations <- 5
parameters <- bednetparams

# Specify a function that runs the run_simulation function for a specified number of iterations:
run_simulations_with_iterations <- function(timesteps, parameters, correlations = NULL, iterations) {

  # Open a list to store the simulation outputs in:
  simulation_outputs <- list()

  # If a correlations object has been provided by the user:
  if(!is.null(correlations)) {

    # Run the simulation for the specified number of iterations:
    for(i in seq(iterations)) {

      # Run the ith simulation with the i-th parameter list:
      simulation_outputs[[i]] <- malariasimulation::run_simulation(timesteps = timesteps,
                                                                   parameters = parameters,
                                                                   correlations = correlations)


      # Append the iteration to the dataframe for identification:
      simulation_outputs[[i]]$iteration <- i

    }

  # If no correlations object has been provided by the user:
  } else {

    # Run the simulation for the specified number of iterations:
    for(i in seq(iterations)) {

      # Run the ith simulation with the i-th parameter list:
      simulation_outputs[[i]] <- malariasimulation::run_simulation(timesteps = timesteps,
                                                                   parameters = parameters)


      # Append the iteration to the dataframe for identification:
      simulation_outputs[[i]]$iteration <- i

    }

  }

  # Create a combined simulation output dataframe
  simulation_outputs <- dplyr::bind_rows(simulation_outputs)

  # Check that the correct number of simulations has been performed:
  if(length(unique(simulation_outputs$iteration)) != iterations) {
    stop("Error: Incorrect number of iterations simulated")
  }

  # Return the combined dataframe of simulation outputs:
  return(simulation_outputs)

}

# Run the simulations using the updated_run_simulation_with_repetitions() function:
outputs_v2 <- run_simulations_with_iterations(
  timesteps = sim_length,
  parameters = bednetparams,
  correlations = correlations,
  iterations = 5)

# Plot the PfPR[2-10] through time for each iteration:
outputs_v2 |>
  mutate(pfpr = n_detect_lm_0_1825 / n_age_0_1825) |>
  ggplot(aes(x = timestep, y = pfpr, colour = as.factor(iteration))) + geom_line() +
  theme_bw() +
  labs(x = "Time (days)", y = "Prevalence", colour = "Iteration") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))

#----- 4) Comparison -------------------------------------------------------------------------------

cowplot::plot_grid(

  # Plot the PfPR[2-10] through time for each iteration:
  outputs_v1 |>
    mutate(pfpr = n_detect_lm_0_1825 / n_age_0_1825) |>
    ggplot(aes(x = timestep, y = pfpr, colour = as.factor(repetition))) + geom_line() +
    theme_bw() +
    labs(x = "Time (days)", y = "Prevalence", colour = "repetition") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_colour_viridis_d() +
    geom_vline(xintercept = bednetstimesteps, linetype = "dotted"),

  # Plot the PfPR[2-10] through time for each iteration:
  outputs_v2 |>
    mutate(pfpr = n_detect_lm_0_1825 / n_age_0_1825) |>
    ggplot(aes(x = timestep, y = pfpr, colour = as.factor(iteration))) + geom_line() +
    theme_bw() +
    labs(x = "Time (days)", y = "Prevalence", colour = "Iteration") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) +
    scale_colour_viridis_d() +
    geom_vline(xintercept = bednetstimesteps, linetype = "dotted"),

  ncol = 2

)


#--------------------------------------------------------------------------------------------------#

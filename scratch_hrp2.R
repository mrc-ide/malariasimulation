# Scratch script for exploring the HRP2 antigen persistence feature.
# Not part of the package - safe to edit/delete freely.
devtools::load_all()
library(ggplot2)

# ---------------------------------------------------------------------------
# 1. Basic setup: a population with treatment enabled, so we get both
#    treated and untreated clinical infections to compare.
# ---------------------------------------------------------------------------
timesteps <- 365 * 20

params <- get_parameters(list(human_population = 5000))
params <- set_equilibrium(params, init_EIR = 10)
params <- set_drugs(params, list(AL_params))
params <- set_clinical_treatment(params, drug = 1, timesteps = 1, coverages = 0.5)

# HRP2 clearance hazard - tweak this and re-run to see the effect on
# n_hrp2_positive (this is a placeholder default, not calibrated).
#
# HRP2 is tracked with a single binary variable (variables$hrp2: 1 =
# currently detectable, 0 = not) rather than a timestamp. It is set to 1 the
# moment a new clinical/asymptomatic infection starts, and stays 1 with NO
# decay for as long as the individual is still infected (state D/A/U/Tr) -
# HRP2 clearance only starts counting down once the underlying infection has
# actually resolved (state == 'S'), at which point a constant per-timestep
# hazard of 1/hrp2_scale applies (the same exponential/memoryless convention
# used for clearing infections elsewhere in the model, e.g. dd/da/du/dt).
# hrp2_scale is the mean number of timesteps between infection clearance and
# HRP2 clearance, and the SAME value is used for everyone - treated
# individuals end up HRP2 positive for less total time only because their
# underlying infection resolves faster (Tr -> S) than an untreated one
# (D -> A -> U -> S), not because of a different HRP2 hazard.
halflife_to_scale <- function(halflife_days) {
  halflife_days / log(2)
}

params <- set_hrp2_parameters(
  params,
  hrp2_asymptomatic_prob = 1,
  hrp2_scale = halflife_to_scale(halflife_days = 60)
)

# Theoretical survival curve implied by the above: the probability that
# someone is still HRP2 positive x timesteps after their infection actually
# cleared (state == 'S'), independent of running a full simulation. This is
# the same for everyone, regardless of whether the resolved infection was
# treated - the treated/untreated difference lives entirely in how long it
# took to reach S in the first place, not in this curve.
hrp2_survival_curve <- data.frame(x = 0:timesteps)
hrp2_survival_curve$survival <- exp(-hrp2_survival_curve$x / params$hrp2_scale)

p0 <- ggplot(hrp2_survival_curve, aes(x = x, y = survival)) +
  geom_line() +
  labs(
    title = "HRP2 persistence survival curve after infection clearance",
    x = "timesteps since infection cleared (state S)",
    y = "P(still HRP2 positive)"
  ) +
  ylim(0, 1) +
  theme_minimal()
print(p0)

# Age-stratified outputs for treated/untreated infections and hrp2 positivity.
# `prevalence` uses the same age breaks as `hrp2` so the two are directly
# comparable (same age groups, same denominator).
params <- set_epi_outputs(
  params,
  hrp2 = c(0, 5 * 365, 100 * 365),
  prevalence = c(0, 5 * 365, 100 * 365)
)

# ---------------------------------------------------------------------------
# 2. Run and pull out the new columns.
# ---------------------------------------------------------------------------
set.seed(1)
sim <- run_simulation(timesteps = timesteps, parameters = params)

hrp2_cols <- grep("hrp2|treated_infections", names(sim), value = TRUE)
print(hrp2_cols)

# ---------------------------------------------------------------------------
# 3. Sanity checks
# ---------------------------------------------------------------------------
stopifnot(all(sim$n_hrp2_positive <= params$human_population))
stopifnot(all(sim$n_treated_infections + sim$n_untreated_infections <= sim$n_infections))
cat("Sanity checks passed\n")

# ---------------------------------------------------------------------------
# 4. Plots
# ---------------------------------------------------------------------------

# New clinical infections per timestep, treated vs untreated
infections_long <- data.frame(
  timestep = rep(sim$timestep, 2),
  count = c(sim$n_treated_infections, sim$n_untreated_infections),
  type = rep(c("treated", "untreated"), each = nrow(sim))
)

p1 <- ggplot(infections_long, aes(x = timestep, y = count, colour = type)) +
  geom_line(alpha = 0.4) +
  geom_smooth(se = FALSE, span = 0.1) +
  labs(title = "New clinical infections per timestep", y = "count") +
  theme_minimal()
print(p1)

# Total HRP2-positive count over time
p2 <- ggplot(sim, aes(x = timestep, y = n_hrp2_positive)) +
  geom_line() +
  labs(title = "Population currently HRP2 antigen positive", y = "n_hrp2_positive") +
  theme_minimal()
print(p2)

# Age-stratified hrp2 positivity (children vs everyone else)
age_long <- data.frame(
  timestep = rep(sim$timestep, 2),
  count = c(sim$n_hrp2_positive_0_1824, sim$n_hrp2_positive_1825_36499),
  age_group = rep(c("under 5", "5+"), each = nrow(sim))
)

p3 <- ggplot(age_long, aes(x = timestep, y = count, colour = age_group)) +
  geom_line() +
  labs(title = "HRP2 positive by age group", y = "count") +
  theme_minimal()
print(p3)

# Prevalence measured by HRP2 positivity vs by light microscopy detection.
# Uses matching age bands (set in set_epi_outputs above) so both are a
# fraction of the same denominator (n_age_<lower>_<upper>).
prevalence_comparison <- data.frame(
  timestep = rep(sim$timestep, 4),
  prevalence = c(
    sim$n_hrp2_positive_0_1824 / sim$n_age_0_1824,
    sim$n_detect_lm_0_1824 / sim$n_age_0_1824,
    sim$n_hrp2_positive_1825_36499 / sim$n_age_1825_36499,
    sim$n_detect_lm_1825_36499 / sim$n_age_1825_36499
  ),
  method = rep(c("HRP2 positive", "Light microscopy"), each = nrow(sim), times = 2),
  age_group = rep(c("under 5", "5+"), each = nrow(sim) * 2)
)

# Facet by method (not age group) with a free y-scale: HRP2 prevalence
# (a few %) and LM prevalence (tens of %) are very different magnitudes, so
# putting both age groups on a shared y-axis per age-group facet squashes
# the smaller metric flat and makes the age groups look identical when
# they aren't. Faceting by method instead makes the under-5 vs 5+ contrast
# visible for each metric on its own scale.
p3b <- ggplot(prevalence_comparison, aes(x = timestep, y = prevalence, colour = method)) +
  geom_line() +
  facet_wrap(~age_group) +
  labs(
    title = "Prevalence: HRP2 positivity vs light microscopy detection",
    y = "prevalence"
  ) +
  theme_minimal()
print(p3b)

cat("\nMean prevalence by age group and method:\n")
print(aggregate(prevalence ~ age_group + method, prevalence_comparison, mean))


##Only run the following if the simulation is low lift
if(!(timesteps > (365 * 2) | params$human_population > 5000)){
  # ---------------------------------------------------------------------------
  # 5. Quick look at how hrp2_scale changes the steady-state HRP2-positive
  #    fraction - useful for building intuition when calibrating.
  # ---------------------------------------------------------------------------
  scales_to_try <- c(5, 20, 60)
  scale_results <- lapply(scales_to_try, function(sc) {
    p <- set_hrp2_parameters(params, hrp2_scale = sc)
    s <- run_simulation(timesteps = timesteps, parameters = p)
    data.frame(timestep = s$timestep, n_hrp2_positive = s$n_hrp2_positive, hrp2_scale = sc)
  })
  scale_df <- do.call(rbind, scale_results)
  scale_df$hrp2_scale <- factor(scale_df$hrp2_scale)

  p4 <- ggplot(scale_df, aes(x = timestep, y = n_hrp2_positive, colour = hrp2_scale)) +
    geom_line() +
    labs(title = "Effect of hrp2_scale on HRP2 positive count", y = "n_hrp2_positive") +
    theme_minimal()
  print(p4)

  # ---------------------------------------------------------------------------
  # 6. Among people who are HRP2 positive, how many are still actively
  #    infected (D/A/U/Tr, no decay yet) vs already cleared and mid-decay
  #    (state S)? run_simulation()/run_resumable_simulation() only return the
  #    rendered per-timestep dataframe, not per-individual variables - to
  #    inspect variables$hrp2/variables$state directly we replicate the same
  #    setup they use internally (see run_resumable_simulation(), R/model.R)
  #    but keep a live reference to `variables`, which the individual package
  #    mutates in place as the simulation runs.
  # ---------------------------------------------------------------------------
  run_simulation_with_variables <- function(timesteps, parameters) {
    correlations <- get_correlation_parameters(parameters)
    variables <- create_variables(parameters)
    events <- create_events(parameters)
    initialise_events(events, variables, parameters)
    renderer <- individual::Render$new(timesteps)
    populate_incidence_rendering_columns(renderer, parameters)
    attach_event_listeners(events, variables, parameters, correlations, renderer)
    vector_models <- parameterise_mosquito_models(parameters, timesteps)
    solvers <- parameterise_solvers(vector_models, parameters)
    lagged_eir <- create_lagged_eir(variables, solvers, parameters)
    lagged_infectivity <- create_lagged_infectivity(variables, parameters)

    individual::simulation_loop(
      processes = create_processes(
        renderer, variables, events, parameters, vector_models, solvers,
        correlations, lagged_eir, lagged_infectivity, timesteps
      ),
      variables = variables,
      events = events,
      timesteps = timesteps
    )

    list(data = renderer$to_dataframe(), variables = variables)
  }

  set.seed(1)
  result <- run_simulation_with_variables(timesteps, params)

  is_positive <- result$variables$hrp2$get_values() == 1
  state_of_positives <- result$variables$state$get_values()[is_positive]

  cat(
    "\nAmong", sum(is_positive), "HRP2-positive individuals at the end of the simulation,",
    "breakdown by current disease state:\n"
  )
  print(table(state_of_positives))
  cat(
    "\n(D/A/U/Tr = still infected, no HRP2 decay yet;",
    "S = infection cleared, currently decaying towards HRP2-negative)\n"
  )
}

# Scratch script for exploring the new HRP2 antigen persistence feature.
# Not part of the package - safe to edit/delete freely.
devtools::load_all()
library(ggplot2)

# ---------------------------------------------------------------------------
# 1. Basic setup: a population with treatment enabled, so we get both
#    treated and untreated clinical infections to compare.
# ---------------------------------------------------------------------------
timesteps <- 365 * 2

params <- get_parameters(list(human_population = 5000))
params <- set_equilibrium(params, init_EIR = 10)
params <- set_drugs(params, list(AL_params))
params <- set_clinical_treatment(params, drug = 1, timesteps = 1, coverages = 0.5)

# HRP2 clearance hazard - tweak these and re-run to see the effect on
# n_hrp2_positive (these are placeholder defaults, not calibrated).
# Treated infections are given their own (faster) clearance rate here via
# hrp2_treated_same_clearance = FALSE; set it back to TRUE (or drop
# hrp2_shape_treated/hrp2_scale_treated) to have everyone clear at the same
# rate instead.
params <- set_hrp2_parameters(
  params,
  hrp2_asymptomatic_prob = 0.7,
  hrp2_shape = 1, hrp2_scale = 500,
  hrp2_treated_same_clearance = FALSE,
  hrp2_shape_treated = 2, hrp2_scale_treated = 10
)

# Theoretical survival curves implied by the above: the probability that an
# infection from x timesteps ago is still HRP2 positive, independent of
# running a full simulation. Shown separately for each group.
hrp2_survival_curve <- data.frame(
  x = rep(0:(6 * max(params$hrp2_scale, params$hrp2_scale_treated)), 2)
)
hrp2_survival_curve$group <- rep(c("untreated/asymptomatic", "treated"), each = nrow(hrp2_survival_curve) / 2)
hrp2_survival_curve$survival <- c(
  weibull_survival(
    hrp2_survival_curve$x[hrp2_survival_curve$group == "untreated/asymptomatic"],
    params$hrp2_shape, params$hrp2_scale
  ),
  weibull_survival(
    hrp2_survival_curve$x[hrp2_survival_curve$group == "treated"],
    params$hrp2_shape_treated, params$hrp2_scale_treated
  )
)

p0 <- ggplot(hrp2_survival_curve, aes(x = x, y = survival, colour = group)) +
  geom_line() +
  labs(
    title = "HRP2 persistence survival curve by group",
    x = "timesteps since infection",
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

# ---------------------------------------------------------------------------
# 5. Quick look at how hrp2_scale changes the steady-state HRP2-positive
#    fraction - useful for building intuition when calibrating.
# ---------------------------------------------------------------------------
scales_to_try <- c(5, 20, 60)
scale_results <- lapply(scales_to_try, function(sc) {
  p <- set_hrp2_parameters(params, hrp2_shape = 2, hrp2_scale = sc)
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
# 6. Quick look at how hrp2_shape changes the steady-state HRP2-positive
#    fraction - useful for building intuition when calibrating.
# ---------------------------------------------------------------------------
shapes_to_try <- c(1, 2, 4, 8)
shape_results <- lapply(shapes_to_try, function(sh) {
  p <- set_hrp2_parameters(params, hrp2_shape = sh, hrp2_scale = params$hrp2_scale)
  s <- run_simulation(timesteps = timesteps, parameters = p)
  data.frame(timestep = s$timestep, n_hrp2_positive = s$n_hrp2_positive, hrp2_shape = sh)
})
shape_df <- do.call(rbind, shape_results)
shape_df$hrp2_shape <- factor(shape_df$hrp2_shape)

p5 <- ggplot(shape_df, aes(x = timestep, y = n_hrp2_positive, colour = hrp2_shape)) +
  geom_line() +
  labs(title = "Effect of hrp2_shape on HRP2 positive count", y = "n_hrp2_positive") +
  theme_minimal()
print(p5)

# ---------------------------------------------------------------------------
# 7. Disaggregated clearance: treated vs untreated/asymptomatic.
#    `params` (set up at the top of the script) already uses disaggregated
#    clearance (hrp2_treated_same_clearance = FALSE). Compare against a
#    "shared clearance rate" baseline, where treated infections clear at the
#    same rate as everyone else, to see the aggregate effect.
# ---------------------------------------------------------------------------
params_shared <- set_hrp2_parameters(params, hrp2_treated_same_clearance = TRUE)

set.seed(1)
sim_shared <- run_simulation(timesteps = timesteps, parameters = params_shared)
# `sim` (from section 2) already used the disaggregated `params`

clearance_comparison <- data.frame(
  timestep = rep(sim$timestep, 2),
  n_hrp2_positive = c(sim$n_hrp2_positive, sim_shared$n_hrp2_positive),
  scenario = rep(c("disaggregated (faster treated clearance)", "shared clearance rate"), each = nrow(sim))
)

p7 <- ggplot(clearance_comparison, aes(x = timestep, y = n_hrp2_positive, colour = scenario)) +
  geom_line() +
  labs(
    title = "Effect of disaggregating treated HRP2 clearance",
    y = "n_hrp2_positive"
  ) +
  theme_minimal()
print(p7)

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
params <- set_hrp2_parameters(params, hrp2_shape = 2, hrp2_scale = 20)

# Age-stratified outputs for treated/untreated infections and hrp2 positivity.
params <- set_epi_outputs(params, hrp2 = c(0, 5 * 365, 100 * 365))

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

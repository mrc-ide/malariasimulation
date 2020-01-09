get_parameters <- function() {
  timestep_to_day <- 1
  list(
    b0    = 0.590076,
    b1    = 0.5,
    ib0   = 43.8787,
    kb    = 2.15506,
    rd    = 1/5,
    ra    = 1/195,
    ru    = 1/110,
    rt    = 1/5,
    ft    = 1/2,
    av1   = .92,
    av2   = .74,
    av3   = .94,
    cd    = 0.068,
    ct    = 0.021896,
    ca    = 1.82425,
    cu    = 0.00062,
    rm    = 1 / (67.6952 * timestep_to_day),
    rb    = 1 / (10 * 365 * timestep_to_day),
    rc    = 1 / (30 * 365 * timestep_to_day),
    ub    = 1 / 7.19919,
    uc    = 1 / 67.6952,
    a0    = 8 * 365 * timestep_to_day,
    theta0  = .0749886,
    theta1  = .0001191,
    ic0     = 18.02366,
    rho     = .85,
    sigma_squared   = 1.67,
    timestep_to_day = 1
  )
}

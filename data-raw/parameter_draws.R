## code to prepare parameter draws: pf
pd_pf <- read.csv("data-raw/parameter_draws_pf.csv")

# Split individual draws_pf into a list of named parameter lists
draws_pf <- 1:max(pd_pf$draw)
parameter_draws_pf <- list()
for(draw in draws_pf){
  pars <- as.list(pd_pf[pd_pf$draw == draw, "value"])
  names(pars) <- pd_pf[pd_pf$draw == draw, "parameter"]
  parameter_draws_pf[[draw]] <- pars
}
names(parameter_draws_pf) <- draws_pf

usethis::use_data(parameter_draws_pf, overwrite = TRUE)

## code to prepare parameter draws: pv
## this dataset does not match the default parameters, but the code should still be functional.
pd_pv <- read.csv("data-raw/parameter_draws_pv.csv", header = T) |> 
  tibble::rownames_to_column(var = "draw") |> 
  dplyr::mutate(phi0  = phi_D_max,
                phi1 = phi_D_min/phi_D_max) |> 
  dplyr::select(-c(loglike, phi_D_min, phi_D_max)) |> 
  tidyr::pivot_longer(cols = !draw, names_to = "parameter", values_to = "value") |> 
  dplyr::mutate(parameter = dplyr::case_when(
    parameter == "u_par" ~ "ua",
    parameter == "phi_LM_max" ~ "philm_max",
    parameter == "phi_LM_min" ~ "philm_min",
    parameter == "A_LM_50pc" ~ "alm50",
    parameter == "K_LM" ~ "klm",
    parameter == "u_clin" ~ "ud",
    parameter == "phi0" ~ "phi0",
    parameter == "phi1" ~ "phi1",
    parameter == "A_D_50pc" ~ "ic0",
    parameter == "K_D" ~ "kc",
    parameter == "A_d_PCR_50pc" ~ "apcr50",
    parameter == "K_d_PCR" ~ "kpcr",
    parameter == "d_PCR_max" ~ "dpcr_max",
    parameter == "d_LM" ~ "da",
    parameter == "P_MI" ~ "pcm",
    parameter == "d_MI" ~ "rm"
  ),
  draw = as.numeric(draw)) |> 
  as.data.frame()

# Split individual draws_pf into a list of named parameter lists
draws_pv <- 1:max(pd_pv$draw)
parameter_draws_pv <- list()
for(draw in draws_pv){
  pars <- as.list(pd_pv[pd_pv$draw == draw, "value"])
  names(pars) <- pd_pv[pd_pv$draw == draw, "parameter"]
  parameter_draws_pv[[draw]] <- pars
}
names(parameter_draws_pv) <- draws_pv

usethis::use_data(parameter_draws_pv, overwrite = TRUE)

## code to prepare parameter draws

#########################################
## Prepare falciparum  parameter draws ##
#########################################
pd <- read.csv("data-raw/parameter_draws.csv")

# Split individual draws into a list of named parameter lists
draws <- 1:max(pd$draw)
parameter_draws_f <- list()
for(draw in draws){
  pars <- as.list(pd[pd$draw == draw, "value"])
  names(pars) <- pd[pd$draw == draw, "parameter"]
  parameter_draws_f[[draw]] <- pars
}
names(parameter_draws_f) <- draws

usethis::use_data(parameter_draws_f, overwrite = TRUE)

###################################
## Prepare vivax parameter draws ##
###################################

#############################
## To address: These vivax parameter draws no not appear to match the default values
#############################

pdv <- read.csv("data-raw/PvMod_posterior_draws.csv") |>
  dplyr::mutate(phi0 = phi_D_max,
                phi1 = phi_D_min/phi_D_max,
                draw = 1:1000) |> 
    dplyr::select(-c(loglike,phi_D_max,phi_D_min)) |> 
  tidyr::pivot_longer(cols = 1:16, names_to = "parameter", values_to = "value") |> 
  dplyr::mutate(parameter = dplyr::case_when(parameter == "phi_LM_max" ~ "philm_max",
                                             parameter == "phi_LM_min" ~ "philm_min",
                                             parameter == "A_LM_50pc" ~ "alm50",
                                             parameter == "K_LM" ~ "klm",
                                             parameter == "u_clin" ~ "uc",
                                             parameter == "A_D_50pc" ~ "ic0",
                                             parameter == "K_D" ~ "kc",
                                             parameter == "A_d_PCR_50pc" ~ "apcr50",
                                             parameter == "K_d_PCR" ~ "kpcr",
                                             parameter == "d_PCR_max" ~ "dpcr_max",
                                             parameter == "d_LM" ~ "da",
                                             parameter == "P_MI" ~ "pcm",
                                             parameter == "d_MI" ~ "rm",
                                             parameter == "sig_het" ~ "sigma_squared",
                                             .default = as.character(parameter))) |> 
  as.data.frame()

# Split individual draws into a list of named parameter lists
draws_v <- 1:max(pdv$draw)
parameter_draws_v <- list()
for(draw in draws_v){
  pars <- as.list(pdv[pdv$draw == draw, "value"])
  names(pars) <- pdv[pdv$draw == draw, "parameter"]
  parameter_draws_v[[draw]] <- pars
}
names(parameter_draws_v) <- draws_v

usethis::use_data(parameter_draws_v, overwrite = TRUE)

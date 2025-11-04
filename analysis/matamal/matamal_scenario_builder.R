#getting parameters for scenario building
df_setting <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_odin_inputs.rds")

df_phi <- read.csv("analysis/matamal/phi_estimates_bijagos.csv", header = TRUE) %>%
  filter(scenario == "Data")

Q0_vals <- c(df_setting$hbi_lower, df_setting$hbi_median, df_setting$hbi_upper)

Q0_vals_in <- c(min(Q0_vals), df_setting$hbi_median[1], max(Q0_vals))
df_phi_in <- c(round(df_phi$lower,3), round(df_phi$mean,3), round(df_phi$upper))
ivm_cov_in <- unique(c(df_setting$lower_ivm_cov, df_setting$mean_ivm_cov, df_setting$upper_ivm_cov))
retention_net_in <- c(12, 21, 38)*30 #AG work https://www.medrxiv.org/content/10.1101/2025.08.27.25334550v1.full.pdf
#DP_cov_in <- unique(c(df_setting$lower_DP_cov, df_setting$mean_DP_cov, df_setting$upper_DP_cov))
DP_cov_in <- unique(df_setting$mean_DP_cov)

df_scenario <- expand.grid(Q0_vals_in, df_phi_in, ivm_cov_in, retention_net_in, DP_cov_in)
names(df_scenario) <- c("Q0", "bites_Bed", "ivm_cov", "retention_net_in", "DP_cov_in")

saveRDS(df_scenario, file = "analysis/matamal/output/df_scenario_sens.rds")

(nrow(df_scenario)*10)/60 #expected run time in hours
#13 hours for intervention arm, control will just be Q0, bites_Bed and net retention (so 9 runs)

#if include DP cov var
(243*10)/60 #40 hours...

##test scenario for the median ivm coverage and Q0

df_test <- df_scenario %>%
  filter(Q0 == Q0_vals_in[2] & ivm_cov == ivm_cov_in[2]) %>%
  mutate(endec_mu = 0.15, 
         wane_endec = 0.01) #update accordingly for the other parameters

saveRDS(df_test, file = "analysis/matamal/output/df_test_scenario.rds")

#then for the endec_mu and wanes that have been mapped across
endec_mapping <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/matamal-runs/malariasim_endec_killing_inputs.rds")
endec_mapping <- endec_mapping %>%
  separate(scenario, into = c("Q0_chr", "cov_chr"), sep = "_Q0_", remove = FALSE) %>%
  mutate(cov_chr = gsub("_cov", "", cov_chr)) %>%
  mutate(Q0= case_when(Q0_chr == "low" ~ Q0_vals_in[1], 
                       Q0_chr == "medium" ~ Q0_vals_in[2], 
                       Q0_chr == "high" ~ Q0_vals_in[3]), 
         ivm_cov = case_when(cov_chr == "low" ~ ivm_cov_in[1], 
                             cov_chr == "medium" ~ ivm_cov_in[2], 
                             cov_chr == "high" ~ ivm_cov_in[3]))

#map carefully
df_scenario_int <- left_join(df_scenario, endec_mapping, by = c("Q0", "ivm_cov"))

df_scenario_int2 <- df_scenario_int %>%
  rename(mu_endec = best_endec_mu, 
         wane_endec = best_wane) %>%
  select(-c(scenario, errors)) %>%
  distinct()
saveRDS(df_scenario_int2, file = "analysis/matamal/output/df_scenario_int.rds")

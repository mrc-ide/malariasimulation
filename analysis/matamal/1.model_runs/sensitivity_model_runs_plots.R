all_results_intervention <- readRDS("analysis/matamal/1.model_runs/sensitivity_output/intervention/intervention_sensitivity_partial_81.rds")
names(all_results_intervention) #can see the unlabelled...oops!

all_results_control <- readRDS("analysis/matamal/1.model_runs/sensitivity_output/control/control_sensitivity_partial_27.rds")
length(all_results_control)#27
df_scenario_control <- readRDS("analysis/matamal/output/df_scenario_sens.rds") %>%
  select(Q0, bites_Bed, retention_net_in, DP_cov_in) %>%
  distinct()

df_scenario_int <- readRDS("analysis/matamal/output/df_scenario_int.rds")
correct_names <- paste0(
  "Q0_", df_scenario_int$Q0,
  "_phi_", df_scenario_int$bites_Bed,
  "_net_", df_scenario_int$retention_net_in,
  "_ivm_cov_", df_scenario_int$ivm_cov,
  "_endec_mu_", df_scenario_int$mu_endec,
  "_wane_endec_", df_scenario_int$wane_endec
)
length(all_results_intervention)
length(correct_names)

# Assign names
names(all_results_intervention) <- correct_names

all_results_control <- lapply(names(all_results_control), function(nm) {
  df <- all_results_control[[nm]]
  df$scenario <- nm
  df
}) |> bind_rows()

all_results_control <- all_results_control %>%
  unnest(cols = everything())

all_results_intervention <- lapply(names(all_results_intervention), function(nm) {
  df <- all_results_intervention[[nm]]
  df$scenario <- nm
  df
}) |> bind_rows()

all_results_intervention <- all_results_intervention %>%
  unnest(cols = everything())

all_results_intervention <- all_results_intervention %>%
  mutate(scenario_short = str_extract(scenario, ".*_ivm_cov_[0-9\\.]+"))



#dates of prevalence surveys etc####
timestep_from_2012 <- function(year, month, day) {
  # Days in each month (non-leap year)
  month_lengths <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  # Years since 2012 (each year = 365 days)
  years_since <- year - 2012
  days_from_years <- years_since * 365
  
  # Days from previous months in the same year
  if (month > 1) {
    days_from_months <- sum(month_lengths[1:(month - 1)])
  } else {
    days_from_months <- 0
  }
  
  # Days within the current month (start at day 1 → offset 0)
  days_from_days <- day - 1
  
  # Total timestep (timestep = 1 corresponds to 1 Jan 2012)
  timestep <- 1 + days_from_years + days_from_months + days_from_days
  return(timestep)
}

#dates of the prevalence surveys
prev_survey_date <- timestep_from_2012(2018, 11, 1) 
prev_survey_date_2021 <- timestep_from_2012(2021, 11, 1)
prev_survey_date_2022 <- timestep_from_2012(2022, 11, 1)


#info on ITN distributions 
distr_campaign <- 6*30 #distribution in June
distr_years <- c(2014, 2017, 2020) #in these years

#prevalence survey information
prev<- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/baseline_prev_malariasim_params.rds")
prev <- prev %>%
  add_row(year = as.factor(2018), mean_prev = 0.148, lower_prev = 14.5, upper_prev = 14.9) %>%
  mutate(mean_prev = mean_prev*100, 
         lower_prev = lower_prev*100, 
         upper_prev = upper_prev*100)

prev_2018 <- prev %>%
  filter(year == 2018) %>%
  mutate(day = prev_survey_date) 

prev_2021 <- prev %>%
  filter(year == 2021) %>%
  mutate(day = prev_survey_date_2021)

prev_2022 <- prev %>%
  filter(year == 2022) %>%
  mutate(day = prev_survey_date_2022)

#cluster-level prevalence survey information 
df_prev_cluster <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/cluster_prev_clean.rds")

df_prev_cluster <- df_prev_cluster %>%
  mutate(day = case_when(year == 2021 ~ prev_2021$day, 
                         year == 2022 ~ prev_2022$day, 
                         TRUE ~ NA_real_))

df_prev_cluster_control <- df_prev_cluster %>%
  filter(arm == 0)

df_prev_cluster_int <- df_prev_cluster %>%
  filter(arm == 1)

#mda campaign arrows 

#dates of the MDAs
DP_july_21 <- timestep_from_2012(2021, 7, 1)
DP_aug_21 <- timestep_from_2012(2021, 8, 1)
DP_sept_21 <- timestep_from_2012(2021, 9, 1)

DP_july_22 <- timestep_from_2012(2022, 7, 1)
DP_aug_22 <- timestep_from_2012(2022, 8, 1)
DP_sept_22 <- timestep_from_2012(2022, 9, 1)

#DP_mda_y_2021 <- c(y_2021+july, y_2021+aug, y_2021+sept)
DP_mda_y_2021 <- c(DP_july_21, DP_aug_21, DP_sept_21)
DP_mda_y_2022 <- c(DP_july_22, DP_aug_22, DP_sept_22)

mda_events <- c(DP_mda_y_2021, DP_mda_y_2022)


mda_campaign_arrows <- data.frame(
  x = ((mda_events-1)/365)+2012,
  xend = ((mda_events-1)/365)+2012,
  y = rep(0.2*100, length(distr_years)),
  yend = rep(0.17*100, length(distr_years))
)

#ID the median scenario
df_setting <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_odin_inputs.rds")

df_phi <- read.csv("analysis/matamal/phi_estimates_bijagos.csv", header = TRUE) %>%
  filter(scenario == "Data")
phi_B <- round(df_phi$mean,3)
Q0 <- df_setting$hbi_median[1]
ivm_cov <- unique(df_setting$mean_ivm_cov)
retention_net <- c(12, 21, 38)*30 #AG work https://www.medrxiv.org/content/10.1101/2025.08.27.25334550v1.full.pdf
retention_net_med <- retention_net[2]
DP_cov <- round(df_setting$mean_DP_cov[1], 3)

target_scenario <- paste0(
  "Q0_", round(Q0, 3),
  "_phi_", round(phi_B, 3),
  "_net_", retention_net_med,
  "_DP_cov_", DP_cov
)


all_results_control <- all_results_control %>%
  mutate(
    target_scenario = scenario == target_scenario,
    scenario_label = ifelse(target_scenario, "Reference scenario", "Other scenarios")
  )

all_results_control %>%
  filter(target_scenario == TRUE)

target_scenario_int <- paste0(
  "Q0_", round(Q0, 3),
  "_phi_", round(phi_B, 3),
  "_net_", retention_net_med,
  #"_DP_cov_", DP_cov,
  "_ivm_cov_", ivm_cov
)


all_results_intervention <- all_results_intervention %>%
  mutate(
    target_scenario = scenario_short == target_scenario_int,
    scenario_label = ifelse(target_scenario, "Reference scenario", "Other scenarios")
  )


saveRDS(all_results_intervention, file = "analysis/matamal/1.model_runs/sensitivity_output/processed_outputs/intervention.rds")
saveRDS(all_results_control, file = "analysis/matamal/1.model_runs/sensitivity_output/processed_outputs/control.rds")


all_results_intervention %>%
  filter(target_scenario == TRUE)

control_plot <- ggplot(all_results_control, aes(x = ((timestep-1)/365)+2012, y = (n_detect_pcr_0_29200/n_age_0_29200)*100, col = as.factor(scenario)))+
  geom_line()+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("All-age qPCR prevalence (%)")+
  xlab("Year")+
  ggtitle("Control Arm")+
  #add ITN distribution campaign arrows
  annotate("segment", x = distr_years[1]+(distr_campaign/365), 
           xend = distr_years[1]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[2]+(distr_campaign/365), 
           xend = distr_years[2]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[3]+(distr_campaign/365), 
           xend = distr_years[3]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100,
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  #add on cluster-level data for each prevalence survey 
  geom_point(data = df_prev_cluster_control, aes(x = (day/365)+2012, y = all_age_prev), col = "slateblue")+
  #2018 prevalence survey info 
  geom_errorbar(data = prev_2018,
                aes(x = ((day-1)/365)+2012, ymin = lower_prev/100, ymax = upper_prev/100),
                color = "red", width = 0.2,
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2018, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  #2021 prevalence survey info 
  geom_errorbar(data = prev_2021,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2021, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  #2022 prevalence survey info 
  geom_errorbar(data = prev_2022,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2022, aes(x = (day/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
 
  #show MDA distributions
  geom_segment(
    data = mda_campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1, col = "orange"
  )+
  #control the axis
  coord_cartesian(xlim = c(2012, 2024 -0.5)) +
  scale_x_continuous(breaks = 2012:2024)


int_plot <- ggplot(all_results_intervention, aes(x = ((timestep-1)/365)+2012, y = (n_detect_pcr_0_29200/n_age_0_29200)*100, col = as.factor(scenario)))+
  geom_line()+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("All-age qPCR prevalence (%)")+
  xlab("Year")+
  ggtitle("Intervention Arm")+
  #add ITN distribution campaign arrows
  annotate("segment", x = distr_years[1]+(distr_campaign/365), 
           xend = distr_years[1]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[2]+(distr_campaign/365), 
           xend = distr_years[2]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[3]+(distr_campaign/365), 
           xend = distr_years[3]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100,
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  #add on cluster-level data for each prevalence survey 
  geom_point(data = df_prev_cluster_int, aes(x = (day/365)+2012, y = all_age_prev), col = "slateblue")+
  #2018 prevalence survey info 
  geom_errorbar(data = prev_2018,
                aes(x = ((day-1)/365)+2012, ymin = lower_prev/100, ymax = upper_prev/100),
                color = "red", width = 0.2,
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2018, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  #2021 prevalence survey info 
  geom_errorbar(data = prev_2021,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2021, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  #2022 prevalence survey info 
  geom_errorbar(data = prev_2022,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2022, aes(x = (day/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  
  #show MDA distributions
  geom_segment(
    data = mda_campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1, col = "orange"
  )+
  #control the axis
  coord_cartesian(xlim = c(2012, 2024 -0.5)) +
  scale_x_continuous(breaks = 2012:2024)

cowplot::plot_grid(control_plot, int_plot)  


#using the highlighter
control_plot_highlight <- ggplot(all_results_control, aes(x = ((timestep-1)/365)+2012, y = (n_detect_pcr_0_29200/n_age_0_29200)*100, col = as.factor(scenario_label)))+
  geom_line()+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("All-age qPCR prevalence (%)")+
  xlab("Year")+
  ggtitle("Control Arm")+
  #add ITN distribution campaign arrows
  annotate("segment", x = distr_years[1]+(distr_campaign/365), 
           xend = distr_years[1]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[2]+(distr_campaign/365), 
           xend = distr_years[2]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[3]+(distr_campaign/365), 
           xend = distr_years[3]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100,
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  #add on cluster-level data for each prevalence survey 
  geom_point(data = df_prev_cluster_control, aes(x = (day/365)+2012, y = all_age_prev), col = "slateblue")+
  #2018 prevalence survey info 
  geom_errorbar(data = prev_2018,
                aes(x = ((day-1)/365)+2012, ymin = lower_prev/100, ymax = upper_prev/100),
                color = "red", width = 0.2,
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2018, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  #2021 prevalence survey info 
  geom_errorbar(data = prev_2021,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2021, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  #2022 prevalence survey info 
  geom_errorbar(data = prev_2022,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2022, aes(x = (day/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  
  #show MDA distributions
  geom_segment(
    data = mda_campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1, col = "orange"
  )+
  #control the axis
  coord_cartesian(xlim = c(2012, 2024 -0.5)) +
  scale_x_continuous(breaks = 2012:2024)+
  scale_colour_manual(values = c("grey", "purple"))

int_plot_highlight <- ggplot(all_results_intervention, aes(x = ((timestep-1)/365)+2012, y = (n_detect_pcr_0_29200/n_age_0_29200)*100, col = as.factor(scenario_label)))+
  geom_line()+
  theme_bw()+
  theme(legend.position = "none")+
  ylab("All-age qPCR prevalence (%)")+
  xlab("Year")+
  ggtitle("Intervention Arm")+
  #add ITN distribution campaign arrows
  annotate("segment", x = distr_years[1]+(distr_campaign/365), 
           xend = distr_years[1]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[2]+(distr_campaign/365), 
           xend = distr_years[2]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100, 
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  annotate("segment", x = distr_years[3]+(distr_campaign/365), 
           xend = distr_years[3]+(distr_campaign/365), y = 0.26*100, yend = 0.21*100,
           arrow = arrow(length = unit(0.25, "cm")), size = 1.5)+
  #add on cluster-level data for each prevalence survey 
  geom_point(data = df_prev_cluster_control, aes(x = (day/365)+2012, y = all_age_prev), col = "slateblue")+
  #2018 prevalence survey info 
  geom_errorbar(data = prev_2018,
                aes(x = ((day-1)/365)+2012, ymin = lower_prev/100, ymax = upper_prev/100),
                color = "red", width = 0.2,
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2018, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  #2021 prevalence survey info 
  geom_errorbar(data = prev_2021,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2021, aes(x = ((day-1)/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  #2022 prevalence survey info 
  geom_errorbar(data = prev_2022,
                aes(x = (day/365)+2012, ymin = lower_prev, ymax = upper_prev),
                color = "red", width = 0.2, 
                inherit.aes = FALSE, size = 1)+
  geom_point(data = prev_2022, aes(x = (day/365)+2012, y = mean_prev), col = "black", shape = 4, 
             inherit.aes = FALSE)+
  
  #show MDA distributions
  geom_segment(
    data = mda_campaign_arrows,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm")),
    size = 1, col = "orange"
  )+
  #control the axis
  coord_cartesian(xlim = c(2012, 2024 -0.5)) +
  scale_x_continuous(breaks = 2012:2024)+
  scale_colour_manual(values = c("grey", "purple"))

cowplot::plot_grid(control_plot_highlight, int_plot_highlight)

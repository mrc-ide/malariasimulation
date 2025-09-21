require(tidyverse)

#plots showing influence of net retention times and timing of MDA in relation to ITN distr on impact

#figure: 3 panels 
#panel A: dynamics plot (prev) with lines:
          #for baseline (no int)
          #high-phiB and long retention
          #high phiB and long retention and endec early 
          #high phiB and long retention and endec late

#panel B: heatmap with 3 different retention levels on x and 2 timing levels on y, show cases averted

#panel C: heatmap with 3 different retention levels on x and 2 different phi-Bed on y, show cases averted

#panelA ####
baseline <- readRDS("analysis/chapter-int-ecology/run_baseline_sim.rds") %>%
  mutate(n_use_net = 0, 
         scenario = "baseline") #put in to ease the rbind
long_bednet <- readRDS("analysis/chapter-int-ecology/run_bednets_long.rds") %>%
  mutate(scenario = "long_bednets")

medium_bednet <- readRDS("analysis/chapter-int-ecology/run_bednets_medium.rds") %>%
  mutate(scenario = "medium_bednets")

short_bednet <- readRDS("analysis/chapter-int-ecology/run_bednets_short.rds") %>%
  mutate(scenario = "short_bednets")

endec_early_long_bednet <- readRDS("analysis/chapter-int-ecology/run_endec_early_long_bednet.rds") %>%
  mutate(scenario = "endec_early_long_bednets")
endec_late_long_bednet <- readRDS("analysis/chapter-int-ecology/run_endec_late_long_bednet.rds")%>%
  mutate(scenario = "endec_late_long_bednets")

dynamics_df <- do.call("rbind", list(baseline, long_bednet, endec_early_long_bednet, 
                                        endec_late_long_bednet))


pals <- c('#1b9e77','#d95f02','#7570b3','#e7298a')

itn_on <- 365*5
sim_length <- 365*15
net_seq <- seq(itn_on, sim_length, by = 3*365)
bednet_timesteps <- net_seq
bednet_start <- bednet_timesteps[1]

dynamics_df$scenario2 <- factor(dynamics_df$scenario, 
                                levels = c("baseline", "long_bednets", "endec_early_long_bednets", 
                                           "endec_late_long_bednets"))



dynamics_plot <- ggplot(dynamics_df, aes(x = (timestep-bednet_start)/365, y = (n_detect_lm_0_1825/n_age_0_1825)*100, 
                        col = as.factor(scenario2)))+
  geom_line(size = 1)+
  ylim(0, 100)+
  ylab("Slide prevalence (%) in children under 5-years-old")+
  theme_bw(base_size = 14)+
  theme(legend.position = c(0.3, 0.2))+
  xlab("Time (years) since first ITN campaign")+
  scale_colour_manual(values = pals, 
                      name = "Scenario", 
                      labels = c("Baseline", "Long ITN retention", "Early ivermectin MDA, long ITN retention", 
                                 "Late ivermectin MDA, long ITN retention"))+
  coord_cartesian(xlim = c(-0.5, 7.5))+
  geom_vline(xintercept = 0, lty = "dashed", size = 1)+ #start time of ITN distributions
  geom_vline(xintercept = (early_t_start-bednet_start)/365, lty = "dashed", col = pals[3], size = 1)+
  geom_vline(xintercept = (early_t_end-bednet_start)/365, lty = "dashed", col = pals[3], size = 1)+
  geom_vline(xintercept = (late_t_start-bednet_start)/365, lty = "dashed", col = pals[4], size = 1)+
  geom_vline(xintercept = (late_t_end-bednet_start)/365, lty = "dashed", col = pals[4], size = 1)


ggplot(dynamics_df, aes(x = (timestep-bednet_start)/365, y = (n_inc_clinical_0_1825/n_age_0_1825)*1000, 
                        col = as.factor(scenario2)))+
  geom_line(size = 1)+
  #ylim(0, 100)+
  ylab("Clinical incidence per 1000 children \n under 5-years-old")+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom")+
  xlab("Time (years)")+
  scale_colour_manual(values = pals, 
                      name = "Scenario", 
                      labels = c("Baseline", "Long ITN retention", "Early ivermectin MDA, long ITN retention", 
                                 "Late ivermectin MDA, long ITN retention"))+
  coord_cartesian(xlim = c(-0.5, 7.5))+
  geom_vline(xintercept = 0, lty = "dashed") #start time of ITN distributions
  
  
#population ITN usage
#1e5 is human pop size


ggplot(long_bednet, aes(x = timestep/365, y = n_use_net/1e5))+
  geom_line()

ggplot(long_bednet, aes(x = timestep/365, y = n_detect_lm_0_1825/n_age_0_1825))+
  geom_line()

odin_mda_early <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_early.rds") %>%
  filter(ref == 2)
odin_mda_medium <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_medium.rds")
odin_mda_late <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/exploring_interactions/MIM_poster/chapter_plots/df_var1_mda_late.rds")


#heatmap
endec_early_medium_bednet <- readRDS("analysis/chapter-int-ecology/run_endec_early_medium_bednet.rds") %>%
  mutate(scenario = "endec_early_medium_bednet")
endec_early_short_bednet <- readRDS("analysis/chapter-int-ecology/run_endec_early_short_bednet.rds") %>%
  mutate(scenario = "endec_early_short_bednet")
endec_late_medium_bednet <- readRDS("analysis/chapter-int-ecology/run_endec_late_medium_bednet.rds") %>%
  mutate(scenario = "endec_late_medium_bednet")
endec_late_short_bednet <- readRDS("analysis/chapter-int-ecology/run_endec_late_medium_bednet.rds") %>%
  mutate(scenario = "endec_late_short_bednet")

dynamics_df2 <- do.call("rbind", list(baseline, long_bednet, endec_early_long_bednet, 
                                      endec_late_long_bednet, 
                                      medium_bednet, endec_early_medium_bednet, 
                                      endec_late_medium_bednet,
                                      short_bednet, endec_early_short_bednet, endec_late_short_bednet))

dynamics_df2 <- dynamics_df2 %>%
  mutate(net_retention = case_when(grepl("long", scenario) ~ "long (5y)",
                                   grepl("medium", scenario) ~"medium (2y)", 
                                   grepl("short", scenario) ~ "short (6m)", 
                                   scenario == "baseline" ~ "no interventions")) %>%
  mutate(mda_timing = case_when(grepl("early", scenario) ~ "ITN + early MDA (6m after ITN campaign)",
                                grepl("late", scenario) ~"ITN + late MDA (2y after ITN campaign)", 
                                TRUE ~ "ITN only"))


#dynamics_df2$scenario2 <- factor(dynamics_df2$scenario, 
#                                levels = c("baseline", "long_bednets", "endec_early_long_bednets", 
#                                           "endec_late_long_bednets",
#                                           "short_bednet", "endec_early_short_bednet", "endec_late_short_bednet"))
#
retention_mda_dynamics <- ggplot(dynamics_df2, aes(x = (timestep-bednet_start)/365, y = (n_detect_lm_0_1825/n_age_0_1825)*100, 
                        col = as.factor(mda_timing)))+
  geom_line(size = 1)+
  facet_wrap(vars(net_retention))+
  ylim(0, 100)+
  ylab("Slide prevalence (%) in children under 5-years-old")+
  theme_bw(base_size = 14)+
  theme(legend.position = "bottom")+
  xlab("Time (years)")+
  #scale_colour_manual(values = pals, 
  #                    name = "Scenario", 
  #                    labels = c("Baseline", "Long ITN retention", "Early ivermectin MDA, long ITN retention", 
  #                               "Late ivermectin MDA, long ITN retention"))+
  coord_cartesian(xlim = c(-0.5, 7.5))+
  geom_vline(xintercept = 0, lty = "dashed") #start time of ITN distributions



ggplot(endec_early_long_bednet, aes(x = timestep, y = n_detect_lm_0_1825/n_age_0_1825))+
  geom_line()+
  geom_line(data = odin_mda_early, aes(x = t, y = slide_prev0to5), col = "red")

#long retention in black, keeps prevalence lower over time
ggplot(endec_early_long_bednet, aes(x = timestep, y = n_detect_lm_0_1825/n_age_0_1825))+
  geom_line()+
  geom_line(data = endec_early_short_bednet, aes(x = timestep, y = n_detect_lm_0_1825/n_age_0_1825), col = "red")



#time periods 
mda_int <- 30
eff_len <- 23
early_IVM <- 180
early_IVM_begin1 <- bednet_timesteps[2]+early_IVM # 6 months into bednet campaign. IVM distr is after the 2nd ITN campaign
early_IVM_start <- c(early_IVM_begin1, early_IVM_begin1+mda_int, early_IVM_begin1+mda_int+mda_int)

early_t_start <- early_IVM_start[1]
early_t_end <- early_IVM_start[1]+180

late_IVM <- 2*365 #2y after
late_IVM_begin1 <- bednet_timesteps[2]+late_IVM #
late_IVM_start <- c(late_IVM_begin1, late_IVM_begin1+mda_int, late_IVM_begin1+mda_int+mda_int)

late_t_start <- late_IVM_start[1]
late_t_end <- late_IVM_start[1]+180

#counterfactual scenarios####

ggplot(baseline, aes(x = timestep, y = n_inc_clinical_0_1825/n_age_0_1825))+
  geom_line()

#perennial so does not matter, can always use this as baseline
baseline_inc_early <- baseline %>%
  #select(timestep, n_inc_clinical_0_1825/n_age_0_1825, n_age_0_1825) %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= early_t_start & timestep <= early_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "baseline")

ggplot(long_bednet, aes(x = timestep, y = n_inc_clinical_0_1825/n_age_0_1825))+
  geom_line()+
  geom_vline(xintercept = early_t_start, col = "red")+
  geom_vline(xintercept = early_t_end, col = "red")

#long retention, early counterfactual
long_ret_early_cf <- long_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= early_t_start & timestep <= early_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "long_ret_early_cf")

#long retention, late counterfactual
long_ret_late_cf <- long_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= late_t_start & timestep <= late_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "long_ret_late_cf")

medium_ret_early_cf <- medium_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= early_t_start & timestep <= early_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "medium_ret_early_cf")

medium_ret_late_cf <- medium_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= late_t_start & timestep <= late_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "medium_ret_late_cf")
  
short_ret_early_cf <- short_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= early_t_start & timestep <= early_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "medium_ret_early_cf")

short_ret_late_cf <- short_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= late_t_start & timestep <= late_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "short_ret_late_cf")

#interventions####
long_ret_early_MDA <- endec_early_long_bednet %>%
  mutate(clin_inc0to5 =(n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= early_t_start & timestep <= early_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "long_ret_early_MDA", 
         eff_baseline = ((baseline_inc_early$tot_cases - tot_cases)/baseline_inc_early$tot_cases)*100, 
         eff_ret_cf = ((long_ret_early_cf$tot_cases - tot_cases)/long_ret_early_cf$tot_cases)*100)

long_ret_late_MDA <- endec_late_long_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= late_t_start & timestep <= late_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "long_ret_late_MDA",
         eff_baseline = ((baseline_inc_early$tot_cases - tot_cases)/baseline_inc_early$tot_cases)*100,
         eff_ret_cf = ((long_ret_late_cf$tot_cases - tot_cases)/long_ret_late_cf$tot_cases)*100)

medium_ret_early_MDA <- endec_early_medium_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= early_t_start & timestep <= early_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "medium_ret_early_MDA",
         eff_baseline = ((baseline_inc_early$tot_cases - tot_cases)/baseline_inc_early$tot_cases)*100, 
         eff_ret_cf = ((medium_ret_early_cf$tot_cases - tot_cases)/medium_ret_early_cf$tot_cases)*100)

medium_ret_late_MDA <- endec_late_medium_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= late_t_start & timestep <= late_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "medium_ret_late_MDA",
         eff_baseline = ((baseline_inc_early$tot_cases - tot_cases)/baseline_inc_early$tot_cases)*100,
         eff_ret_cf = ((medium_ret_late_cf$tot_cases - tot_cases)/medium_ret_late_cf$tot_cases)*100)

short_net_early_MDA <- endec_early_short_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= early_t_start & timestep <= early_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "short_net_early_MDA", 
         eff_baseline = ((baseline_inc_early$tot_cases - tot_cases)/baseline_inc_early$tot_cases)*100,
         eff_ret_cf = ((short_ret_early_cf$tot_cases - tot_cases)/short_ret_early_cf$tot_cases)*100)

short_net_late_MDA <- endec_late_short_bednet %>%
  mutate(clin_inc0to5 = (n_inc_clinical_0_1825/n_age_0_1825)*1000) %>%
  filter(timestep >= late_t_start & timestep <= late_t_end) %>%
  summarise(tot_cases = sum(clin_inc0to5)) %>%
  mutate(scenario = "short_net_late_MDA", 
         eff_baseline = ((baseline_inc_early$tot_cases - tot_cases)/baseline_inc_early$tot_cases)*100,
         eff_ret_cf = ((short_ret_late_cf$tot_cases - tot_cases)/short_ret_late_cf$tot_cases)*100)

net_mda_df <- do.call("rbind", list(long_ret_early_MDA, long_ret_late_MDA, 
                                    medium_ret_early_MDA, medium_ret_late_MDA, 
                                    short_net_early_MDA, short_net_late_MDA))

net_mda_df <- net_mda_df %>%
  mutate(net_retention = case_when(grepl("long", scenario) ~ "long (5y)",
                           grepl("medium", scenario) ~"medium (2y)", 
                           grepl("short", scenario) ~ "short (6m)")) %>%
  mutate(mda_timing = case_when(grepl("early", scenario) ~ "early (6m)",
                                grepl("late", scenario) ~"late (2y)"))


heatmap_mda_duration <- ggplot(net_mda_df, aes(x = mda_timing, y = net_retention, fill = eff_baseline))+
  geom_tile()+
  geom_text(aes(label = paste0(round(eff_baseline, 1), "%")),
            col = "white", size = 5)+
  theme_bw(base_size = 14)+
  scale_fill_viridis_c(name = "Cases averted (%) in \n under 5-year-olds due to endectocide
                       from start of MDA to 6m later")+
  xlab("Timing of MDA in relation to ITN campaign")+
  ylab("Net retention period")+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

heatmap_mda_duration_eff_ret_cf <- ggplot(net_mda_df, aes(x = mda_timing, y = net_retention, fill = eff_ret_cf))+
  geom_tile()+
  geom_text(aes(label = paste0(round(eff_ret_cf, 1), "%")),
            col = "white", size = 5)+
  theme_bw(base_size = 14)+
  scale_fill_viridis_c(name = "Cases averted (%) in \n under 5-year-olds due to endectocide
                       from start of MDA to 6m later")+
  xlab("Timing of MDA in relation to most recent ITN campaign")+
  ylab("Net retention period")+
  theme(legend.position = "bottom", 
        legend.direction = "horizontal")

net_retention_figure <- cowplot::plot_grid(dynamics_plot, heatmap_mda_duration_eff_ret_cf, 
                   labels = c("A", "B"))
ggsave(net_retention_figure, file = "analysis/chapter-int-ecology/plots/net_retention_figure_3.pdf")

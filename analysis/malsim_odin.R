#comparing runs from odin against malariasim

hannah_int <- readRDS("W:/endectocides-cluster/raw_outputs/odin-fits/output_mod1.rds")
nilani_int <- readRDS("W:/endectocides-cluster/raw_outputs/odin-fits/odin_exp_decay_compare.rds")

hannah_no_int <- readRDS("W:/endectocides-cluster/raw_outputs/odin-fits/output_mod1_no_int.rds")
nilani_no_int <- readRDS("W:/endectocides-cluster/raw_outputs/odin-fits/odin_exp_decay_compare_no_int.rds")

#no other interventions####
#what is the best-fitting endec_mu and wane_endec for init_EIR = 100?

time_1 <- nilani_no_int %>%
  filter(t == 1) 
time_1$EIR_tot

devtools::load_all()
year <- 365
sim_length <- year*10
human_population <- 100
starting_EIR <- 120

simparams <- get_parameters(overrides = list(
  human_population = human_population,
  endec = TRUE,
  prevalence_rendering_min_ages = 0,    
  prevalence_rendering_max_ages = 5 * 365, 
  individual_mosquitoes = FALSE))


simparams <- set_equilibrium(parameters = simparams, init_EIR = starting_EIR) ##is this the right order? Should we set to eqm later on?
mosq_params <- gamb_params
mosq_params$Q0 <- 0.94
simparams <- set_species(simparams, species = list(mosq_params), 
                         proportions = c(1))

IVM_begin1 <- year*6
mda_int <- 30
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)
eff_len <- 23
steps <- (IVM_start[1]):(IVM_start[3]+eff_len)
#steps <- (IVM_start[1]): (IVM_start[2])
#coverage_in <- c(rep(0.9, eff_len), rep(0, mda_int-eff_len), 
#                 rep(0.9, eff_len), rep(0, mda_int-eff_len), 
#                 rep(0.9, eff_len+1)) #90% coverage for first 23 days of each month


#endec_on <- c(rep(IVM_start[1], eff_len), rep(0, mda_int - eff_len),
#            rep(IVM_start[2], eff_len), rep(0, mda_int - eff_len), 
#            rep(IVM_start[3], eff_len+1))

endec_on <- c(as.numeric(IVM_start[1]):as.numeric(IVM_start[1]+eff_len), 
              as.numeric(IVM_start[2]):as.numeric(IVM_start[2]+eff_len), 
              as.numeric(IVM_start[3]):as.numeric(IVM_start[3]+eff_len))
endec_ts <- c(IVM_start[1], IVM_start[2], IVM_start[3])

endec_params <- set_endectocide(parameters = simparams, timesteps = steps,
                                endec_on = endec_on, endec_ts = endec_ts) 

#output_no_int <- run_simulation_with_repetitions(timesteps = sim_length, overrides =  endec_params, repetitions = 5)
output_no_int <- run_simulation(timesteps = sim_length, parameters = endec_params)



write_rds(output_no_int, file = "analysis/output_endec_no_int.rds")

output_endec <- readRDS("analysis/output_endec_no_int.rds")

output_endec_summary <- output_endec %>%
  group_by(timestep) %>%
  summarise(mean_mu_gamb = mean(mu_gamb), 
           mean_M_gamb = mean(total_M_gamb), 
           mean_EIR_gamb = mean(EIR_gamb), 
           mean_prev0to5 = mean(n_detect_lm_0_1825/n_age_0_1825))

ggplot(output_endec_summary, aes(x = timestep/365, y = mean_prev0to5))+
  geom_line()
 

ggplot(output_endec_summary, aes(x = timestep/365, y = mean_mu_gamb))+
  geom_line()+
  xlim(5, 7.5)

range(output_endec_summary$mean_mu_gamb)

#read in the odin models
odin_models <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/exploring_interactions/malariasim-odin/odin_compare_no_int.rds")
covs <- unique(odin_models$ivm_cov) #0 is for endec_mu model 

#do for high coverage first
odin_compare_highcov <- odin_models %>%
  filter(ivm_cov == covs[2])



odin_models %>%
  filter(model == "odin_exp_decay") %>%
  ggplot()+
  aes(x = t, y = slide_prev0to5)+
  geom_line()+
  facet_wrap(vars(init_EIR, Q0))

ggplot(odin_compare_highcov, aes(x = t, y = slide_prev0to5, col = as.factor(model)))+
  geom_line()+
  facet_wrap(vars(init_EIR, Q0))+
  ylim(0, 1)

output_endec2 <- output_endec_summary %>%
  select(timestep, mean_M_gamb, mean_EIR_gamb, mean_prev0to5) %>%
  mutate(model_type = "malariasim") %>%
  rename(t = timestep, 
         mv = mean_M_gamb, 
         EIR_tot = mean_EIR_gamb, 
         slide_prev0to5 = mean_prev0to5)



output_endec3 <- output_endec2 %>%
  rename(slide_prev0to5_malsim = slide_prev0to5) %>%
  select(t, slide_prev0to5_malsim)

ggplot(odin_models, aes(x = t, y = slide_prev0to5, col = as.factor(model)))+
  geom_point()

odin_models_wide <- odin_models %>%
  select(t, slide_prev0to5, model) %>%
  pivot_wider(names_from = model, values_from = slide_prev0to5) %>%
  rename(slide_prev0to5_hazards = hazards, 
         slide_prev0to5_endec_mu_decay = odin_exp_decay)

output_join <- left_join(output_endec3, odin_models3, by = "t") #join the malariasim and odin models

#compare the model predictions in the time of endectocide killing
output_join_compare <- output_join %>%
  filter(t >= IVM_start[1] & t <= IVM_start[3] + 23) %>%
  select(t, slide_prev0to5_hazards, slide_prev0to5_malsim)

lm_compare <- lm(slide_prev0to5_hazards ~ slide_prev0to5_malsim, data = output_join_compare)
summary(lm_compare) # P value <2e-16 adj Rsq 0.9998

confint(lm_compare, level = 0.95)

ggplot(output_join_compare, aes(x = slide_prev0to5_malsim, y = slide_prev0to5_hazards))+
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) + # y = x line
  theme_bw() +
  labs(
    title = "Comparison of Slide Prevalence: Malsim vs Odin",
    x = "Slide Prevalence U5s (Malsim)",
    y = "Slide Prevalence U5s (Odin hazards)"
  ) +
  annotate("text", x = min(output_join_compare$slide_prev0to5_malsim), 
           y = max(output_join_compare$slide_prev0to5_hazards) + 0.02, 
           label = paste0("Adj. R² = ", signif(summary(lm_compare)$adj.r.squared, digits = 4)), 
           hjust = 0, vjust = 1, color = "blue")+
  xlim(0.55, 0.85)+
  ylim(0.55, 0.85)


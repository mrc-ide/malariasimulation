#first comparing the effects of adding the waning effect to mortality in malariasim vs the other models

devtools::load_all()

require(tidyverse)
year <- 365
sim_length <- 3 * year
human_population <- 1e5
starting_EIR <- 120

require(foresite)

UGA <- foresite::UGA$seasonality

UGA <- UGA %>%
  filter(name_1 == "Tororo")


#simparams <- get_parameters(overrides = list(
#  human_population = human_population,
#  endec = TRUE,
#  prevalence_rendering_min_ages = 0,    
#  g0 = UGA$g0,
#  g = c(UGA$g1,UGA$g2,UGA$g3),
#  h = c(UGA$h1, UGA$h2, UGA$h3),
#  prevalence_rendering_max_ages = 5 * 365, 
#  individual_mosquitoes = FALSE))

simparams <- get_parameters(overrides = list(
  human_population = human_population,
  endec = TRUE,
  prevalence_rendering_min_ages = 0,    
  prevalence_rendering_max_ages = 5 * 365, 
  individual_mosquitoes = FALSE))


simparams <- set_equilibrium(parameters = simparams, init_EIR = starting_EIR) ##is this the right order? Should we set to eqm later on?
mosq_params <- gamb_params
mosq_params$Q0 <- 0.9
simparams <- set_species(simparams, species = list(mosq_params), 
                         proportions = c(1))

IVM_begin1 <- 365+100 # d100 in
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
                                 endec_on = endec_on, endec_ts = endec_ts) #Error in if (expected_bites > 0) { : 


#bednetstimesteps <- c(1, 4) * year # The bed nets will be distributed at the end of the first and the 4th year. 
#
#bednetparams <- set_bednets(
#  endec_params,
#  timesteps = bednetstimesteps,
#  coverages = c(.5, .5),  # Each round is distributed to 50% of the population.
#  retention = 5 * year, # Nets are kept on average 5 years
#  dn0 = matrix(c(.533, .533), nrow = 2, ncol = 1), # Matrix of death probabilities for each mosquito species over time
#  rn = matrix(c(.56, .56), nrow = 2, ncol = 1), # Matrix of repelling probabilities for each mosquito species over time
#  rnm = matrix(c(.24, .24), nrow = 2, ncol = 1), # Matrix of minimum repelling probabilities for each mosquito species over time
#  gamman = rep(2.64 * 365, 2) # Vector of bed net half-lives for each distribution timestep
#)


#run simulation with no interventions 

output_no_int <- run_simulation(timesteps = sim_length, parameters = simparams)

output_no_int %>%
  filter(timestep == 400) %>%
  summarise(prev = n_detect_lm_0_1825/n_age_0_1825)

#output_no_int <- run_simulation_with_repetitions(timesteps = sim_length, overrides = simparams, repetitions = 30)
#write_rds(output_no_int, file = "analysis/output_no_int")

output_endec <- run_simulation_with_repetitions(timesteps = sim_length, overrides = endec_params, repetitions = 30)
write_rds(output_endec, file = "analysis/output_endec.rds")


output_endec_summary <- output_endec %>%
  group_by(timestep) %>%
  summarise(mean_mu_gamb = mean(mu_gamb), 
            mean_M_gamb = mean(total_M_gamb), 
            mean_EIR_gamb = mean(EIR_gamb), 
            mean_prev0to5 = mean(n_detect_lm_0_1825/n_age_0_1825))
output_endec_summary$mean_prev0to5



ggplot(output_endec, aes(x = timestep, y = mu_gamb, col = as.factor(repetition)))+
  geom_line()+
 # geom_hline(yintercept = 0.132, col = "red")+
  geom_vline(xintercept = IVM_start[1], col = "red", linetype = "dashed")+
  geom_vline(xintercept = IVM_start[2], col = "red", linetype = "dashed")+
  geom_vline(xintercept = IVM_start[3], col = "red", linetype = "dashed")+
  theme_bw()

ggplot(output_endec, aes(x = timestep, y = total_M_gamb, col = as.factor(repetition)))+
  geom_line()+
  geom_vline(xintercept = IVM_start[1], col = "red")+
  geom_vline(xintercept = IVM_start[2], col = "red")+
  geom_vline(xintercept = IVM_start[3], col = "red")

ggplot(output_endec, aes(x = timestep, y = EIR_gamb, col = as.factor(repetition)))+
  geom_line()+
  geom_hline(yintercept = 0.132, col = "red")

ggplot(output_endec, aes(x = timestep, y = n_detect_lm_0_1825/n_age_0_1825, 
                         col = as.factor(repetition)))+
  geom_line()+
  ylim(0, 0.8) #malariasim slightly overestimates. 

#read in the odin models
odin_models <- readRDS("C:/Users/nc1115/Documents/github/ivRmectin/analysis/exploring_interactions/malariasim-odin/hazards_endec_mu_compare.rds")
covs <- unique(odin_models$ivm_cov)
odin_models <- odin_models %>%
  filter(ivm_cov == covs[2]) %>%
  select(t, mv, EIR_tot, slide_prev0to5, model_type)

ggplot(odin_models, aes(x = t, y = slide_prev0to5, col = as.factor(model_type)))+
  geom_line()

output_endec2 <- output_endec_summary %>%
  select(timestep, mean_M_gamb, mean_EIR_gamb, mean_prev0to5) %>%
  mutate(model_type = "malariasim") %>%
  rename(t = timestep, 
         mv = mean_M_gamb, 
         EIR_tot = mean_EIR_gamb, 
         slide_prev0to5 = mean_prev0to5)

head(output_endec2)
head(odin_models)

output <- rbind(output_endec2, odin_models)

plot(output_endec2$slide_prev0to5, odin_models$slide_prev0to5) #do model stats to calculate R squared. 

output_endec3 <- output_endec2 %>%
  rename(slide_prev0to5_malsim = slide_prev0to5) %>%
  select(t, slide_prev0to5_malsim)

odin_models3 <- odin_models %>%
  select(t, slide_prev0to5, model_type) %>%
  pivot_wider(names_from = model_type, values_from = slide_prev0to5) %>%
  rename(slide_prev0to5_hazards = hazards, 
         slide_prev0to5_endec_mu_decay = endec_mu_decay)

output_join <- left_join(output_endec3, odin_models3, by = "t") #join the malariasim and odin models


#compare the model predictions in the time of endectocide killing
output_join_compare <- output_join %>%
  filter(t >= IVM_start[1] & t <= IVM_start[3] + 23) %>%
  select(t, slide_prev0to5_hazards, slide_prev0to5_malsim)

lm_compare <- lm(slide_prev0to5_hazards ~ slide_prev0to5_malsim, data = output_join_compare)
summary(lm_compare) # P value <2e-16 adj Rsq 0.9998

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

prev_plot <- ggplot(output, aes(x = t, y = slide_prev0to5, col = as.factor(model_type)))+
  geom_line(linewidth = 1, alpha = 0.5)+
  theme_bw()+
  ylim(0, 1)

init_prev <- output %>%
  filter(t == 400) %>%
  group_by(model_type) %>%
  summarise(init_prev = slide_prev0to5)



out_summary <- output %>%
  filter(t >= IVM_start[1] & t <= IVM_start[3] + 23) %>%
  group_by(model_type) %>%
  summarise(mean_prev = mean(slide_prev0to5, na.rm = TRUE)) 


out_summary_join <- left_join(out_summary, init_prev)

ggplot(out_summary_o)

out_summary_join %>%
  mutate(abs_red_prev = init_prev - mean_prev, 
         rel_red_prev = ((init_prev - mean_prev)/init_prev)*100)

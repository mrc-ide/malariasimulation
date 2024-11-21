devtools::load_all()

require(tidyverse)
year <- 365
sim_length <- 8 * year
human_population <- 10000
starting_EIR <- 65 #this seems to hit the correct initial prevalence
mu_endec_high_cov <- 0.16 #when cov is 90%
mu_endec_low_cov <- 0.06

simparams <- get_parameters(overrides = list(
  human_population = human_population,
  endec = TRUE,
  prevalence_rendering_min_ages = 0,    
  mu_endec = 0.16,
  prevalence_rendering_max_ages = 5 * 365, 
  individual_mosquitoes = FALSE))



#steps <- (3*365):(4*365)
#make this coverage 23 days on and off from 8.5y into simulation 


IVM_begin1 <- net_seq[3]+180 # 6 months into new net distribution
mda_int <- 30
IVM_start <- c(IVM_begin1, IVM_begin1+mda_int, IVM_begin1+mda_int+mda_int)
eff_len <- 23
steps <- (IVM_start[1]):(IVM_start[3]+eff_len)
#steps <- (IVM_start[1]): (IVM_start[2])
#steps <- (22*365-44):(25*365-44)


coverage_in <- c(rep(0.9, eff_len), rep(0, mda_int-eff_len), 
                 rep(0.9, eff_len), rep(0, mda_int-eff_len), 
                 rep(0.9, eff_len+1)) #90% coverage for first 23 days of each month
#coverage_in <- c(rep(0.9, eff_len), rep(0, mda_int - eff_len+1))
#coverage_in <- c(rep(1, length(steps)))
#mosq_params <- gamb_params
#mosq_params$Q0 <- 0.92

#mosq_params_zoo <- gamb_params
#mosq_params_zoo$Q0 <- 0.2

#mosq_params$Q0 <- gamb_params$Q0*1/3 #when add this in, the mu goes really low.....

#simparams <- set_species(simparams, species = list(mosq_params), 
 #                        proportions = c(1))


endec_params <- set_endectocide(parameters = simparams, timesteps = steps, coverages =coverage_in) #Error in if (expected_bites > 0) { : 

simparams <- set_equilibrium(parameters = endec_params, init_EIR = starting_EIR) ##is this the right order? Should we set to eqm later on?


output <- run_simulation(timesteps = sim_length, parameters = simparams)


ggplot(output, aes(x = timestep, y = mu_gamb))+
  geom_line()

ggplot(output, aes(x = timestep, y = total_M_gamb))+
  geom_line()+
  ylim(0, 50e4)

#get a drop in mv when we have Q0*cov, but when we do Q0*1/3*cov, get an increase in mv.Even if we artificially reduce Q0 by a factor of 1/3, doesn't work
#think we don't need to account for the 1 bite every 3 days because it is captured through tau1 and tau2
ggplot(output, aes(x = timestep, y = n_detect_0_1825 /n_0_1825))+
  geom_line()+ 
  ylim(0.2, 1)

#malariasim slightly overestimates the drop in prevalence. 

ggplot(output, aes(x = timestep, y = mu_gamb))+
  geom_line()

#out_odin <- readRDS(file = "../ivRmectin/analysis/exploring_interactions/hazards_mu_h_combined.rds")
#out_Q0 <- unique(out_odin$Q0)
#out_cov <- unique(out_odin$ivm_cov)

odin_mods <- readRDS(file = "../ivRmectin/analysis/exploring_interactions/malariasim-odin/mod_compare.rds")
odin_mods <- odin_mods %>%
  filter(model_type != "mu_h") #weird model, take it out
out_cov <- unique(odin_mods$ivm_cov)
out_odin_test <- odin_mods %>%
  filter(ivm_cov == out_cov[2])

#the endec_mu model slightly overestimates the impact of ivermectin
ggplot(out_odin_test, aes(x = t, y= mv, col = as.factor(model_type)))+
  geom_line()

ggplot(out_odin_test, aes(x = t, y = EIR_tot, col = as.factor(model_type)))+
  geom_line()

out_odin_test %>%
  filter(model_type == "hazards") %>%
  ggplot(aes(x = t, y = EIR_tot, col = as.factor(model_type)))+
  geom_line()+
  ylim(0, 105)

output$model_type <- "malariasim"

output_test <- output %>%
  select(timestep, EIR_gamb, n_detect_0_1825, n_0_1825, model_type, total_M_gamb) %>%
  mutate(Q0 = 0.9, ivm_cov = 0.9) %>%
  rename(t = timestep, 
         EIR_tot = EIR_gamb, 
         mv = total_M_gamb) %>%
  mutate(slide_prev0to5 = n_detect_0_1825/n_0_1825) %>% 
  select(-c(n_detect_0_1825, n_0_1825))

out_odin_test <- out_odin_test %>%
  select(t, EIR_tot, model_type, Q0, ivm_cov, slide_prev0to5, mv)


comp_models <- rbind(output_test, out_odin_test)

ggplot(output_test, aes(x = t, y = mv, col = as.factor(model_type)))+
  geom_line()

ggplot(comp_models, aes(x = t, y = slide_prev0to5, col = as.factor(model_type)))+
  geom_line()+
  ylim(0, 1)+
  theme_minimal()+
  labs(col = "Model type")

ggplot(comp_models, aes(x = t, y = mv, col = as.factor(model_type)))+
  geom_line()+
  theme_minimal()+
  labs(col = "Model type")

ggplot(comp_models, aes(x = t, y = slide_prev0to5, col = as.factor(model_type)))+
  geom_line(size = 1)+
  theme_minimal()+
  ylab("Slide prevalence under 5s")+
  labs(col = "Model type")+
  ylim(0, 0.75)+
  scale_color_manual(values = c("red", "green", "blue"), 
                     labels = c("Odin version of malariasim (model C)", 
                                "Hannah & Charlie's model (model A)", 
                                "Malariasimulation with endectocide (model D)"))+
  theme(legend.position = c(0.4, 0.5))


ggplot(comp_models, aes(x = t, y = EIR_tot, col = as.factor(model_type)))+
  geom_line(size = 1)+
  theme_minimal()+
  ylab("EIR")+
  labs(col = "Model type")


#is the absolute drop in prevalence same for the different models?

write_rds(comp_models, file = "../ivRmectin/analysis/exploring_interactions/malariasim-odin/comp_models.rds")

comp_models %>%
  dplyr::filter(between(t, 1000,2000)) %>%
  select(t, model_type, slide_prev0to5) %>%
  spread(key = model_type, value = slide_prev0to5) %>%
  summarise(mean_prev_malariasim = mean(malariasim), 
            mean_prev_hazards = mean(hazards), 
            mean_prev_endec_mu = mean(endec_mu))


comp_models %>%
  filter(between(t, IVM_start1[1], IVM_start3[3]+23)) %>%
  spread(key = model_type, value = slide_prev0to5) %>%
  summarise(mean_prev_malariasim = mean(malariasim), 
            mean_prev_hazards = mean(hazards), 
            mean_prev_endec_mu = mean(endec_mu)) %>%
  mutate(rel_red_prev = m)
  
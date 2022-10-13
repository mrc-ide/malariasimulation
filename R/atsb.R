# install.packages("~/Documents/GitHub/malariasimulation_Ace_params-master/",
#                  repos=NULL,
#                  type="source")

setwd("~/Documents/Imperial/GitHub")
devtools::install_github("pahanc/malariasimulation_Ace_params")
options(repos = c(
  mrcide = 'https://mrc-ide.r-universe.dev',
  CRAN = 'https://cloud.r-project.org'))

# Install some packages
install.packages("malariasimulation")
library(malariasimulation)

setwd("GitHub/malariasimulation/")
library(devtools)
load_all(".")

devtools::install_github("nmoghaddas19/malariasimulation")

11# environment(custom_death_rate) <- asNamespace('malariasimulation')
# assignInNamespace("death_rate", custom_death_rate, ns = "malariasimulation")

Sites_fitted_2021_03_09 <- readRDS("~/GitHub/Sites_fitted_2021_03_09.RDS")
pyre_params <- read.csv("~/Github/pyrethroid_only_nets.csv")
mali <- Sites_fitted_2021_03_09[which(Sites_fitted_2021_03_09$NAME_0=="Mali"),]
rownames(mali) <- 1:nrow(mali)
zambia <- Sites_fitted_2021_03_09[which(Sites_fitted_2021_03_09$NAME_0=="Zambia"),]
rownames(zambia) <- 1:nrow(zambia)
kenya <- Sites_fitted_2021_03_09[which(Sites_fitted_2021_03_09$NAME_0=="Kenya"),]
rownames(kenya) <- 1:nrow(kenya)

year <- 365
month <- 30
sim_length <- 5 * year
human_population <- 1000
starting_EIR <- 50

# Seasonality ----
mali_vectors <- matrix(0, nrow = ncol(mali), ncol = length(mali$vectors[[1]]))
mali_bednets <- matrix(0, nrow = ncol(mali), ncol = length(mali$interventions[[1]]))
mali_seasonality
mali_resistance
for (i in 1:nrow(mali)){
 mali_vectors[,i] <- mali$vectors
}

mali_params <- get_parameters(
  list(
    human_population = human_population,
    model_seasonality = TRUE, # Let's try a bi-modal model
    g0 = 0.28605,
    g = c(0.20636, -0.0740318, -0.0009293),
    h = c(0.173743, -0.0730962, -0.116019),
    mu_atsb = 0.10
  )
)

simparams3 <- set_bednets(simparams3, c(365,730.5), c(0.5,0.5),
                          retention = 5 * year,
                          dn0 = matrix(c(.533, .45), nrow=2, ncol=1),
                          rn = matrix(c(.56, .5), nrow=2, ncol=1),
                          rnm = matrix(c(.24, .24), nrow=2, ncol=1),
                          gamman = rep(2.64 * 365, 2))

atsb_params <- set_atsb(mali_params, 365:1095, rep(0.9,731))
# atsb_params <- set_atsb(simparams3, c(365,730), c(0.8,0.7))

p <- set_equilibrium(atsb_params, starting_EIR)

output_atsb <- run_simulation(sim_length, p)
plot(1:sim_length, output$n_detect_730_3650/output$n_730_3650, pch=20, cex=0.7,
     ylim=c(0,1), frame.plot = F, xlab="Time", ylab="Prevalence", col="black")
points(1:sim_length, output_atsb$n_detect_730_3650/output_atsb$n_730_3650, pch=20, cex=0.7,
       ylim=c(0,1),  xlab="Time", ylab="Prevalence", col="red")
abline(v=c(1,3)*year, col="red")
#plot(1:sim_length, output$n_infections)
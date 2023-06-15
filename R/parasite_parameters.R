#' @title Parameterise the parasite species to use in the model
#'
#' @description
#' The `get_parameters` function assumes as default that Plasmodium falciparum is being modelled (parasite = "falciparum").
#' `set_parasite` takes the P. falciparum parameter set and alters it to reflect vivax parameters when the parasite = "vivax" argument is given.
#' If parasite = "falciparum" is given, it simply returns the given parameter set.
#' The remaining documentation describes the changes that get made when "vivax" is selected.
#' 
#' fixed state transition changes:
#' * dt - the delay for humans to move from state Tr to Ph; default = 1 (falc. = 5)
#' * da - the delay for humans to move from state A to U; default = 10 (falc. = 195)
#' * du - the delay for humans to move from state U to S (removed: now modelled by an equation, see "duration of subpatent infection")
#' 
#' duration of subpatent infection (new parameters to calculate "du"):
#' * du_max - Maximum duration of subpatent infection: default = 70 (new)
#' * du_min - Minimum duration of subpatent infection: default = 10 (new)
#' * ku - Shape parameter: default = 4.602 (new)
#' * au50 - Scale parameter: default = 9.9 (new)
#' 
#' immunity decay rates changes
#' * rm - decay rate for maternal immunity to clinical disease; default = 35.148 (falc. = 67.6952)
#' 
#' immunity decay rates changes (pre-erythrocitic/blood and severe immunity no longer modelled)
#' * rvm - decay rate for maternal immunity to severe disease (removed)
#' * rva - decay rate for acquired immunity to severe disease (removed)
#' * rb - decay rate for acquired pre-erythrocytic (removed)
#' 
#' probability of pre-erythrocytic infection (now a parameter: b)
#' * b - probability of pre-erythrocytic infection: default = 0.5 (new: replaces equation)
#' * b0  - maximum probability due to no immunity (removed)
#' * b1 - maximum reduction due to immunity (removed)
#' * ib0 - scale parameter (removed)
#' * kb - shape parameter (removed)
#' 
#' probability of clinical infection changes
#' * phi0 - maximum probability due to no immunity; default = 0.8605 (falc. = 0.791666)
#' * phi1 - maximum reduction due to immunity; default = 0.018 (falc. = 0.000737)
#' * ic0 - scale parameter; default = 11.538 (falc. = 18.02366)
#' * kc - shape parameter; default = 2.25 (falc. = 2.36949)
#' 
#' probability of severe infection (no longer modelled)
#' * theta0 - maximum probability due to no immunity (removed)
#' * theta1 - maximum reduction due to immunity (removed)
#' * kv - scale parameter (removed)
#' * fv0 - shape parameter (removed)
#' * av - age dependent modifier (removed)
#' * gammav - age dependent modifier (removed)
#' * iv0 - age dependent modifier (removed)
#'
#' immunity reducing probability of detection (differently modelled)
#' * fd0 - time-scale at which immunity changes with age (removed)
#' * ad - scale parameter relating age to immunity (removed)
#' * gammad - shape parameter relating age to immunity (removed)
#' * d0 - maximum probability due to immunity; default = 0.8918 (new)
#' * d1 - minimum probability divided by maximum due to immunity; default = 0.0043/0.8918 (falc. = 0.160527)
#' * id0 - scale parameter; default = 27.52 (falc. = 1.577533)
#' * kd - shape parameter; default = 2.403 (falc. = 0.476614)
#' 
#' immunity boost grace period changes
#' * uc - period in which clinical immunity cannot be boosted; default = 7.85 (falc. = 6.06349)
#' * ud - period in which immunity to detectability cannot be boosted; default = 19.77 (falc. = 9.44512)
#' 
#' immunity boost grace periods (pre-erythrocytic/blood and severe immunity not modelled)
#' * ub - period in which pre-erythrocytic immunity cannot be boosted (removed)
#' * uv - period in which severe immunity cannot be boosted (removed)
#' 
#' infectivity towards mosquitoes changes
#' * cd - infectivity of clinically diseased humans towards mosquitoes; default = 0.8 (falc. = 0.068)
#' * cu - infectivity of sub-patent infection; default = 0.035 (falc. = 0.0062)
#' * ct - infectivity of treated infection; default = 0.4 (falc. = 0.022)
#' * gamma1 - parameter for infectivity of asymptomatic humans (removed: asymptomatic infectivity now a parameter: see "ca")
#' * ca - infectivity of asymptomatic infection; default = 0.1 (new: replaces age-dependent parameterisation)
#' 
#' unique biting rate changes
#' * sigma_squared - heterogeneity parameter; default = 1.29 (falc. = 1.67)
#' 
#' incubation period changed
#' * de - Duration of the human latent period of infection; default = 10 (falc. = 12)
#' * delay_gam - Lag from parasites to infectious gametocytes; default = 0 (falc. = 12.5)
#' * dem - Extrinsic incubation period in mosquito population model; default = 8.4 (falc. = 10)
#' 
#' human mortality parameter changes
#' * pcm - new-born clinical immunity relative to mother's; default = 0.421 (falc. = 0.774368)
#' * pvm - new-born severe immunity relative to mother's (removed: severe immunity not modelled)
#' 
#' initial immunity values changes
#' * init_ivm - the immunity from severe disease at birth (removed: severe immunity no longer modelled)
#' * init_ib - the initial pre-erythrocitic immunity (removed: pre-erythrocytic/blood immunity no longer modelled)
#' * init_iva - the initial acquired immunity from severe disease (removed: severe immunity no longer modelled)
#' * init_idm - the initial maternal immunity to detectability: default = 0 (new: now modelled)
#'  
#'  hypnozoite parameters introduced
#' * f - relapse rate: default = 0.024 (new)
#' * gammal - clearance rate: default = 0.0026 (new)
#'  
#'  miscellaneous parameter changes
#' * parasite - Plasmodium species; default = "falciparum"
#' 
#' 
#' @param parameters the model parameters
#' @param parasite parasite species name: "falciparum" or "vivax"
#' @export
set_parasite <- function(parameters, parasite) {
  if (length(parasite) != 1) {
    stop('Only one parasite species can be selected')
  }
  if (!parasite %in% c("vivax","falciparum")){
    stop('Parasite species must be "vivax" or "falciparum"')
  }
  if (parasite == "falciparum"){
    return(parameters)
  }
  
  if(parasite == "vivax"){
    
    # Under misc
    parameters$parasite <- "vivax"
    
    ## fixed state transitions
    # dd is the same
    parameters$dt <- 1
    parameters$da <- 10
    # du is now an equation
    # parameters <- parameters[names(parameters) != "du"]
    # del, dl, dpl, mup and mum are the same
    
    ## Duration of subpatent infection
    parameters$du_max <- 70 # Maximum duration of subpatent infection
    parameters$du_min <- 10 # Minimum duration of subpatent infection
    parameters$ku <- 4.602  # Shape parameter
    parameters$au50 <- 9.9 # Scale parameter
      
    ## Immunity decay rates
    parameters$rm <- 35.148
    # rvm, rva: severe parameters and immunity do not exist in the vivax model
    # rb: blood (pre-erythrocytic) immunity does not exist in the vivax model
    # parameters <- parameters[!names(parameters) %in% c("rvm","rb","rva")]
    # rc and rid are the same
    
    ## Probability of pre-erythrocytic infection
    # Now represented by a parameter (b)
    # parameters <- parameters[!names(parameters) %in% c("b0","b1","ib0","kb")]
    parameters$b <- 0.5

    ## Probability of clinical infection
    parameters$phi0 <- 0.8605  # Equivalent to the maximum probability of clinical infection
    parameters$phi1 <- 0.018/0.8605 # Equivalent to the minimum probability of clinical infection divided by the maximum
    parameters$ic0 <- 11.538 # Scale parameter
    parameters$kc <- 2.25 # Shape parameter

    ## Probability of severe infection
    # Severe infection is not modelled in vivax
    # parameters <- parameters[!names(parameters) %in% c("theta0","theta1","kv","fv0","av","gammav","iv0")]
    
    ## Immunity reducing probability of detection
    # This function is no longer age dependent
    # parameters <- parameters[!names(parameters) %in% c("fd0","ad","gammad")]
    parameters$d0 <- 0.8918 # Equivalent to maximum probability of detection
    parameters$d1 <- 0.0043/0.8918 # Equivalent to minimum probability of detection divided by the maximum
    parameters$id0 <- 27.52 # Scale parameter
    parameters$kd <- 2.403# Shape parameter
    
    ## Immunity boost grace periods
    parameters$uc <- 7.85 # period in which clinical immunity cannot be boosted
    parameters$ud <- 19.77 # period in which immunity to detectability cannot be boosted
    # Immunity to pre-erythrocity and severe infection are not modelled in vivax
    # parameters <- parameters[!names(parameters) %in% c("ub","uv")]
    
    ## Infectivity towards mosquitoes
    parameters$cd <- 0.8 # diseased
    parameters$cu <- 0.035 # sub-patent infection
    parameters$ct <- 0.4 # treated
    # Infectivity of asymptomatic disease is now a parameter: ca, not based on gamma1
    # parameters <- parameters[!names(parameters) %in% c("gamma1")]
    parameters$ca <- 0.1 # asymptomatic

    ## Unique biting rate
    # a0 and rho are the same
    parameters$sigma_squared <- 1.29
    
    ## Incubation period
    parameters$de <- 10 # latent period
    parameters$delay_gam <- 0 # lag from parasites to infectious gametocytes
    parameters$dem <- 8.4 # extrinsic incubation period
    
    # Human mortality parameters
    # average_age is the same
    parameters$pcm <- 0.421
    # maternal severe immunity (pvm) is not modelled
    # parameters <- parameters[!names(parameters) %in% c("pvm")]
    
    ## Seasonality and carrying capacity parameters
    # These parameters (model_seasonality, g0, g, h, gamma, rainfall_floor) are mosquito associated and remain the same.
    
    ## Larval mortality
    # These parameters (me, ml) are mosquito associated and remain the same.
    
    ## Initial state proportions
    # Remain the same
    
    ## Initial immunity values
    # Blood and severe immunity removed
    # parameters <- parameters[!names(parameters) %in% c("init_ivm","init_ib","iva")]
    # Maternal immunity to detectable disease is modelled in vivax
    parameters$init_idm <- 0

    ## Vector biology
    # These parameters (beta, total_M, init_foim, species, species_proportions, blood_meal_rates, Q0, foraging time) are mosquito associated and remain the same.

    ## Feeding cycle
    # These parameters (bednets, phi_bednets, k0, spraying, phi_indoors) are mosquito associated and remain the same.
    
    ## Treatment parameters
    # These parameters are switched off as default

    ## PEV, MDA, SMC, PMC, TBV parameters
    # Switched off as default
    
    ## Rendering 
    # These parameters stay the same
    # Although I could consider adding hypnozoite number within these... Or average hypnozoite number...

    ## Miscellaneous 
    # These parameters stay the same
    
    ## Hypnozoite parameters
    parameters$f <- 0.024 # relapse rate
    parameters$gammaL <- 0.0026 # clearance rate

  }
  
  parameters
}

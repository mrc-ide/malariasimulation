# malariasimulation 1.6.0

  * Fix MDA bug where undetectable asymptomatics are treated
  * New vignettes
  * Progress bar for long simulations
  * Individual mosquitoes off by default
  * New vaccination code:
    * pre-erythrocytic vaccine functions have been renamed to pev
    * pev functions have PEVProfiles for alternate pev vaccines and boosters
  * Specify carrying capacity over time

# malariasimulation 1.4.0

  * Treatment number rendering for all treatments not just effective
  * Default rendering for drug-based interventions
  * Bed nets and IRS continue post individual's death
  * ITN scheduling does not overwrite on receipt of >1 net

# malariasimulation 1.3.0

  * Custom demography
  * Rainfall floor
  * Documentation fixes for vector controls
  * Simple model for severe disease
  * Alignment of immunity function with legacy model

# malariasimulation 1.2.0

  * added a `news.md` file to track changes to the package.
  * n_inc_clinical includes treated humans
  * disaggregate EIR and state counts by mosquito species
  * remove redundant Total_M and and global EIR outputs
  * new vignette on how to calibrate EIR to prevalence

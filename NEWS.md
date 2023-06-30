# malariasimulation 1.6.0

  * New vignettes
  * Stephensi parameters
  * Fix MDA bug where undetectable asymptomatics are treated

# malariasimulation 1.5.0

  * New vaccination code:
    * pre-erythrocytic vaccine functions have been renamed to pev
    * pev functions have PEVProfiles for alternate pev vaccines and boosters

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

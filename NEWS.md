# malariasimulation 2.0.0

Release v2.0.0. This will mark our shift to [semver 2.0.0](https://semver.org/spec/v2.0.0.html) for versioning, i.e:

> Given a version number MAJOR.MINOR.PATCH, increment the:
> 
>  1. MAJOR version when you make incompatible API changes
>  2. MINOR version when you add functionality in a backward compatible manner
>  3. PATCH version when you make backward compatible bug fixes
> 
> Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

## Transparent changes

* Fix progress bar bug by @pwinskill in https://github.com/mrc-ide/malariasimulation/pull/258
* Fix bug with pev min_wait by @giovannic in https://github.com/mrc-ide/malariasimulation/pull/264
* Fix competing hazards for mass and EPI pev: by @giovannic in https://github.com/mrc-ide/malariasimulation/pull/265
* Improve correlation tests. by @plietar in https://github.com/mrc-ide/malariasimulation/pull/287
* Assign names to all processes. by @plietar in https://github.com/mrc-ide/malariasimulation/pull/325
* Correcting default average_age parameter by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/326
* Add a range check to bitset_index. by @plietar in https://github.com/mrc-ide/malariasimulation/pull/331
* Rewrite the exponential decay process in C++. by @plietar in https://github.com/mrc-ide/malariasimulation/pull/285
* recovery_rates variable has been renamed as (disease) progression_rates. by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/334

### Transparent changes to model behaviour

* Moved mortality to the end of processes to resolve order of competing… by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/292
* calculate_eir function correction by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/296
* Older age groups (older than equilibrium age object) are not correctl… by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/301
* Competing hazards recovery infection resolution by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/300
* Bug/maternal immunity sampling by @tbreweric in https://github.com/mrc-ide/malariasimulation/pull/279

### Improved documentation

* Update description of the dt parameter. by @plietar in https://github.com/mrc-ide/malariasimulation/pull/288
* Clarify some intervention parametrization documentation. by @plietar in https://github.com/mrc-ide/malariasimulation/pull/291
* Parameter and output documentation. by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/299
* Hex by @pwinskill in https://github.com/mrc-ide/malariasimulation/pull/336

## New features

 * Pausing and resuming of simulations by @plietar in https://github.com/mrc-ide/malariasimulation/pull/280 https://github.com/mrc-ide/malariasimulation/pull/286 https://github.com/mrc-ide/malariasimulation/pull/293 https://github.com/mrc-ide/malariasimulation/pull/323
* Function for setting age-structured outputs. by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/297
* Add GitHub workflow for continuous benchmarking. by @plietar in https://github.com/mrc-ide/malariasimulation/pull/304
* adding in n_spray output by @htopazian in https://github.com/mrc-ide/malariasimulation/pull/310
* Antimalarial resistance by @tbreweric in https://github.com/mrc-ide/malariasimulation/pull/263 https://github.com/mrc-ide/malariasimulation/pull/282 https://github.com/mrc-ide/malariasimulation/pull/289
* Add R21 PEV profiles by @lmhaile in https://github.com/mrc-ide/malariasimulation/pull/333

## Breaking changes

* Custom carrying capacity is set as a relative scaler, not absolute value by @pwinskill in https://github.com/mrc-ide/malariasimulation/pull/272
* Implement time varying coverage for PEV boosters: by @giovannic in https://github.com/mrc-ide/malariasimulation/pull/270
* Add asymmetric transmission in mixing in metapopulation modelling: by @giovannic in https://github.com/mrc-ide/malariasimulation/pull/277 https://github.com/mrc-ide/malariasimulation/pull/259
* lm pcr prevalence explicit by @RJSheppard in https://github.com/mrc-ide/malariasimulation/pull/298

## New Contributors
* @lmhaile made their first contribution in https://github.com/mrc-ide/malariasimulation/pull/333

**Full Changelog**: https://github.com/mrc-ide/malariasimulation/compare/v1.6.0...v2.0.0

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

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mda_parameters.R
\name{set_smc}
\alias{set_smc}
\title{Parameterise a Seasonal Malaria Chemoprevention}
\usage{
set_smc(parameters, drug, start, end, frequency, min_age, max_age, coverage)
}
\arguments{
\item{parameters}{a list of parameters to modify}

\item{drug}{the index of the drug to administer}

\item{start}{the timestep to start}

\item{end}{the last timestep for the intervention}

\item{frequency}{the number of timsteps between doses}

\item{min_age}{the minimum age of the target population exclusive (in timesteps)}

\item{max_age}{the maximum age of the target population exclusive (in timesteps)}

\item{coverage}{the proportion of the target population who will recieve the
drug}
}
\description{
Parameterise a Seasonal Malaria Chemoprevention
}
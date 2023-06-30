[![R build status](https://github.com/mrc-ide/malariasimulation/workflows/R-CMD-check/badge.svg)](https://github.com/mrc-ide/malariasimulation/actions)
[![codecov](https://codecov.io/github/mrc-ide/malariasimulation/branch/master/graphs/badge.svg)](https://codecov.io/github/mrc-ide/malariasimulation)

# malariasimulation

Imperial College London's next malaria simulation. The main goals are make an
extensible, maintainable and fast simulation to evaluate and report on malaria
intervention strategies.

The model is defined in this package, whereas the simulation is executed using
the [individual](https://github.com/mrc-ide/individual) package.

## Installation

Please note, malariasimulation is only compatible with R >= 4.0.0

You can install the binary package on Windows and Mac using the following
command:

```R
# Enable this universe
options(repos = c(
    mrcide = 'https://mrc-ide.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

# Install some packages
install.packages('malariasimulation')
```

The package can also be installed from github using the "remotes" library. Note, this
method requires [RBuildTools](https://cran.r-project.org/bin/windows/Rtools/)

```R
library('remotes')
install_github('mrc-ide/malariasimulation')
```

For development it is most convenient to run the code from source. You can
install the dependencies in RStudio by opening the project and selecting "Build" >
"Install and Restart"

Command line users can execute:

```R
library('remotes')
install_deps('.', dependencies = TRUE)
```

Docker users can achieve the same with one line

```bash
docker build . -f docker/Dockerfile.dev -t [your image name]
```

## Usage

To run the malaria simulation with the default parameters for 100 timesteps, you
can execute the following code:

```R
library('malariasimulation')

output <- run_simulation(100)
```

Please see [vignettes](https://mrc-ide.github.io/malariasimulation/articles/Model.html) for more detailed use.

## Code organisation

*model.R* - is the entry point for the model. It creates the different
components and passes them to the `individual` package for simulation.

*variables.R* - contains the specific variables we track for each individual.

*processes.R* - defines the changes in individuals over time. It collects
process functions from "infection.R", "mosquito_emergence.R" and "mortality.R"

*tests* - are divided into unit and integration tests. Integration tests are
strongly recommended for large process functions and unit tests for model
calculations.

## Contributing

Pull requests are welcome. For major changes, please open an issue first to
discuss what you would like to change.

Please make sure to update tests and documentation as appropriate.

Code reviews will be carried out in-line with RESIDE-IC's [review
process](https://reside-ic.github.io/articles/pull-requests/)

## License
[MIT](https://choosealicense.com/licenses/mit/)

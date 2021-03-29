# covidestim <img src="man/figures/logo.png" width="120" align="right" />

Real-time estimates of the true size and trajectory of local COVID-19 epidemics
are key metrics to guide policy responses. *covidestim* is a Bayesian
nowcasting approach that explicitly accounts for reporting delays and secular
changes in case ascertainment to generate real-time estimates of COVID-19
epidemiology on the basis of reported cases and deaths. Using this approach, we
can estimate time trends in infections, symptomatic cases, and deaths for all
50 US states and the District of Columbia from early-March through the present.
The *covidestim* package contains all you need to reproduce our work, but it's
really targeted at users who want to try running our model on their own data.
As the package allows for different kinds of input data, custom priors, and two
different fitting algorithms, *covidestim* is appropriate for:

- Efficient real-time estimation of key timeseries at the U.S. state- and county-level
  (these results are available on [covidestim.org](https://covidestim.org))

- Retrospective modeling of a specific geographic area for which you have
  special insight. For example, *covidestim* can help you model a county where
  you happen to have line-list case and death data by date of occurrence, as
  well as clinical data informing your belief about the rates at which patients
  transition from symptomatic to severe (hospitalizable), and from severe to
  dead.

- Exploratory modeling of non-US geographies. *covidestim*'s extensive
  customizability allows you to leverage beliefs about delay structures, case
  ascertainment, and more in order to set up a model tailored to your
  epidemiologic setting.

*covidestim* provides a well-documented interface for all phases of the modeling
process: adding input data, specifying priors, running the model, producing 
summary timeseries and graphs.

# Installation

```r
# Currently, Covidestim is not availble on CRAN
# install.packages("covidestim")

# Install development version from GitHub. This requires that the 'devtools'
# package be installed.
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("covidestim/covidestim")
```

`rstan`, *covidestim*'s key dependency, has a dependency on Google's V8 JS
library. V8 often proves troublesome to install. If you experience problems,
see [this GitHub issue](https://github.com/stan-dev/rstan/issues/831) for
solutions. Or, just use our Docker container, linked below.

# Usage

```r
library(covidestim)

cfg <- covidestim(ndays = 120, region = 'New York', pop_size = get_pop('New York')) +
  input_cases( example_nyc_data('cases')) +
  input_deaths(example_nyc_data('deaths'))

result <- run(cfg)

resultSummary <- summary(result)
resultGraphs  <- viz(result)
```

# Resources

* [Read our preprint](https://www.medrxiv.org/content/10.1101/2020.06.17.20133983v1)
* [Ask a question](mailto:marcus.russi@yale.edu?subject=covidestim)
* [Use our Docker image](https://hub.docker.com/repository/registry-1.docker.io/covidestim/covidestim), an easy way to get started without dealing with dependencies, and an easy way to use *covidestim* in an HPC or cloud environment
* [Open an issue](https://github.com/covidestim/covidestim/issues) (GitHub
  issues for bug reports, feature requests)

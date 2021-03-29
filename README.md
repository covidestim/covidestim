# covidestim <img src="man/figures/logo.png" width="120" align="right" />

Real-time estimates of the true size and trajectory of local COVID-19 epidemics
are key metrics to guide policy responses. *covidestim* is a Bayesian
nowcasting approach that explicitly accounts for reporting delays and secular
changes in case ascertainment to generate real-time estimates of COVID-19
epidemiology on the basis of reported cases and deaths. Using this approach, we
can estimate time trends for a number of important outcomes:

- Infections (including infections that were never diagnosed)
- Symptomatic cases 
- Deaths
- Effective reproduction number (R<sub>t</sub>)
- ...and more

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

- Exploratory modeling of non-U.S. geographies. *covidestim*'s extensive
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

# Install the development version from GitHub. This requires that the
# 'devtools' package be installed.
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

# Set up a model on 120 days of example NYC data
cfg <- covidestim(ndays = 120, region = 'New York', pop_size = get_pop('New York')) +
  input_cases(example_nyc_data('cases')) +
  input_deaths(example_nyc_data('deaths'))

result <- run(cfg)

# Get a data.frame of all timeseries outcomes, and some overview plots
resultSummary <- summary(result)
resultGraphs  <- viz(result)
```

# Changelog

_covidestim_ has been under active development since March 2020. In that time,
we have made large structural changes to the core model, added and removed data
sources, tweaked priors, added a new fitting method, and more. We are
constantly working to provide an implementation that is consistent with
the best available evidence on the epidemiology of COVID-19. To that end, the
`master` branch generally sees a merge every 1-3 months. While we do not
currently use a versioning system to catalogue these changes, we maintain a
detailed history of major changes [here](https://www.covidestim.org/updates.pdf).

# Resources

* _[covidestim/covidestim-sources](https://github.com/covidestim/covidestim-sources)_
  provides utilities to easily clean up-to-date case/death data for U.S. states
  and counties, plus Puerto Rico and D.C.
* _[covidestim/dailyFlow](https://github.com/covidestim/dailyFlow)_ is a
  [Nextflow](https://nextflow.io/) routine for batch execution of _covidestim_
  in an HPC or cloud environment.
* [Read our preprint](https://www.medrxiv.org/content/10.1101/2020.06.17.20133983v1)
* [Ask a question](mailto:marcus.russi@yale.edu?subject=covidestim)
* [Use our U.S. state/county model results](https://covidestim.org), updated
  daily. The latest results can be downloaded as a `.csv` from a stable URL.
* [Use our Docker image](https://hub.docker.com/repository/registry-1.docker.io/covidestim/covidestim),
  an easy way to get started without dealing with dependencies, and an easy way
  to use *covidestim* in an HPC or cloud environment.
* [Open an issue](https://github.com/covidestim/covidestim/issues) (GitHub
  issues for bug reports, feature requests)

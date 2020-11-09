# covidestim <img src="man/figures/logo.png" width="120" align="right" />

Real-time estimates of the true size and trajectory of local COVID-19 epidemics are
key metrics to guide policy responses. *covidestim* is a Bayesian nowcasting approach that
explicitly accounts for reporting delays and secular changes in case ascertainment to generate
real-time estimates of COVID-19 epidemiology on the basis of reported cases and deaths. Using
this approach, we can estimate time trends in infections, symptomatic cases, and deaths for all 50 US
states and the District of Columbia from early-March through the present. The *covidestim* package
contains all you need to reproduce our work, but it's really targeted at users who want to try
running our model on their own state- or county-level data.

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

# Usage

```r
library(covidestim)

cfg <- covidestim(ndays=80) + input_cases(df) + input_deaths(df2, type = "occurred")

result <- run(cfg)

result.summary <- summary(result)
result.graphs  <- viz(result)
```

# Resources

* [Read our preprint](https://www.medrxiv.org/content/10.1101/2020.06.17.20133983v1)
* [Ask a question](mailto:marcus.russi@yale.edu?subject=covidestim)
* [Use our Docker image](https://hub.docker.com/repository/registry-1.docker.io/covidestim/covidestim)
* [Open an issue](https://github.com/covidestim/covidestim/issues) (GitHub
  issues for bug reports, feature requests)

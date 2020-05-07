# Covidcast

Quick description of Covidcast and what it does

# Installation

```r
# Currently, Covidcast is not availble on CRAN
# install.packages("covidcast")

# Install development version from GitHub. This requires that the 'devtools'
# package be installed.
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("mel-hc/covidcast")
```

# Usage

- To see how to run `covidcast` on example NYC data, see `vignette("nyc")`

- To see how to run `covidcast` on your own data, see `vignette("ILINet")`.

# Resources

- [Ask a question](mailto:marcus.russi@yale.edu?subject=covidcast)
- [Open an issue](https://github.com/mel-hc/covidcast/issues) (GitHub
  issues for bug reports, feature requests)

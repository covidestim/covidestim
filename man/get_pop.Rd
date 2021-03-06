% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covidestim.R
\name{get_pop}
\alias{get_pop}
\title{Population estimates for US states and counties}
\usage{
get_pop(region)
}
\arguments{
\item{region}{A string with the state name, or the FIPS code}
}
\value{
State/county population as a numeric, or an error
}
\description{
Returns 2019 census estimate of state or county population
}
\examples{
get_pop('Connecticut')
get_pop('09009')

}

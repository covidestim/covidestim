#!/usr/bin/env Rscript

# This file needs to be run from the package root directory!

source('data-raw/generate_nyc_daily_counts.R')
source('data-raw/generate_popsize.R')
source('data-raw/generate_ifrs.R')

ifr_state  <- ifrs$state
ifr_county <- ifrs$county
ifr_prior  <- ifrs$ifr_prior

usethis::use_data(
  ifr_state, ifr_county, ifr_prior, nyc_data, pop_state, pop_county,
  internal = TRUE, overwrite = TRUE
)

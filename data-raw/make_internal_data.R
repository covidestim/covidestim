#!/usr/bin/env Rscript

# This file needs to be run from the package root directory!

source('data-raw/generate_nyc_daily_counts.R')
source('data-raw/generate_popsize.R')
source('data-raw/generate_ifrs.R')
source('data-raw/generate_logor_vac.R')

ifr_state  <- ifrs$state
ifr_county <- ifrs$county

usethis::use_data(
  ifr_state, ifr_county, nyc_data, pop_state, pop_county,
  logor_vac_state, logor_vac_county,
  internal = TRUE, overwrite = TRUE
)

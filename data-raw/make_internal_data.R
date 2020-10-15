#!/usr/bin/Rscript

# This file needs to be run from the package root directory!

source('data-raw/generate_nyc_daily_counts.R')
source('data-raw/generate_popsize.R')
source('data-raw/generate_ifrs.R')

ifr_state  <- ifrs$ifr_state
ifr_county <- ifrs$ifr_county

usethis::use_data(
  ifr_state, ifr_county, nyc_data, pop_state, pop_county,
  internal = TRUE, overwrite = TRUE
)

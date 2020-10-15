#!/usr/bin/Rscript

source('data-raw/generate_nyc_daily_counts.R')
source('data-raw/generate_ifrs.R')

ifr_state  <- ifrs$ifr_state
ifr_county <- ifrs$ifr_county

usethis::use_data(
  ifr_state, ifr_county, nyc_data,
  internal = TRUE, overwrite = TRUE
)

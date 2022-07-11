#!/usr/bin/env Rscript

# This file needs to be run from the package root directory!

source('data-raw/generate_ct_weekly_counts.R')
source('data-raw/generate_popsize.R')
source('data-raw/generate_init_immunity.R')
source('data-raw/generate_ifrs.R')

ifr_state  <- ifrs$state
ifr_county <- ifrs$county

usethis::use_data(
  ifr_state, ifr_county, ct_data, pop_state, pop_county,
  imm_county, imm_state,
  internal = TRUE, overwrite = TRUE
)

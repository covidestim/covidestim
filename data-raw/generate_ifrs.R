#!/usr/bin/Rscript
library(readr)

ifrs <- local({

  ifr_county <- read_csv(
    'data-raw/ifr-data/ifr_county_12-2-2020.csv',
    col_types = 'cccnn'
  )

  ifr_state <- read_csv(
    'data-raw/ifr-data/ifr_state_date_OR_12-2-2020.csv',
    col_types = cols(.default = col_number(), Date = col_date())
  )

  list(
    county = ifr_county,
    state  = ifr_state
  )
});

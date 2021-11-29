#!/usr/bin/Rscript
library(readr)

ifrs <- local({

  ifr_county <- read_csv(
    'data-raw/ifr-data/ifr_county_11-17-2021.csv',
    col_types = 'cccn'
  )

  ifr_state <- read_csv(
    'data-raw/ifr-data/ifr_state_12-12-2020.csv',
    col_types = 'cnn'
  )
  
  list(
    county    = ifr_county,
    state     = ifr_state
  )
});

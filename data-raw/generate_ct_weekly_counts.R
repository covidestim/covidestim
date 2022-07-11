#!/usr/bin/Rscript

ct_data <- local({
  library(magrittr)
  
  readr::read_csv(
    'data-raw/ct_weekly_counts.csv',
    col_types = readr:: cols(
      date  = readr:: col_date(format = "%Y-%m-%d"),
      .default = readr::col_number()
    )
  ) %>% 
    tidyr::replace_na(list(cases=0, deaths=0, RR=1,boost=0,hosp=0))
})

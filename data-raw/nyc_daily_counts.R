readr::read_csv(
  'data-raw/nyc_daily_counts.csv',
  col_types = cols(date     = readr::col_date(format = "%m/%d/%y"),
                   .default = readr::col_number())
) %>%
tidyr::replace_na(
  list(cases=0, hosps=0, deaths=0)
) -> nyc_data

usethis::use_data(nyc_data, internal = TRUE)

library(magrittr)

readr::read_csv(
  'data-raw/nyc_daily_counts.csv',
  col_types = readr:: cols(
    date     = readr::col_date(format = "%m/%d/%Y"),
    .default = readr::col_number()
  )
) %>%
tidyr::replace_na(
  list(cases=0, deaths=0)
) -> nyc_data

usethis::use_data(nyc_data, internal = TRUE, overwrite = TRUE)

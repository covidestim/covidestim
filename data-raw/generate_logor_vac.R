logor_vac_county <- readr::read_csv(
  'data-raw/vac-data/logor-vac-fips.csv',
  col_type = 'cnn'
)

logor_vac_state <- readr::read_csv(
  'data-raw/vac-data/logor-vac-state.csv',
  col_type = 'cnn'
)

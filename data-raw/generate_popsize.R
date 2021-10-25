pop_state  <- readr::read_csv(
  'data-raw/pop-data/statepop.csv',
  col_types = 'cn'
)

pop_county <- readr::read_csv(
  'data-raw/pop-data/fipspop.csv',
  col_types = 'cn'
)

pop_state_age <- readr::read_csv(
  'data-raw/pop-data/statepop-age.csv',
  col_types = 'cn'
)

pop_county_age <- readr::read_csv(
  'data-raw/pop-data/fipspop-age.csv',
  col_types = 'cn'
)

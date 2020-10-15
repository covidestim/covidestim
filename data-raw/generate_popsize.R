pop_state  <- readr::read_csv(
  'data-raw/pop-data/statepop.csv',
  col_types = 'cn'
)

pop_county <- readr::read_csv(
  'data-raw/pop-data/fipspop.csv',
  col_types = 'cn'
)

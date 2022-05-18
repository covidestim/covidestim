imm_state  <- readr::read_csv(
  'data-raw/imm-data/state-imm-init.csv',
  col_types = 'cnn'
)

imm_county <- readr::read_csv(
  'data-raw/imm-data/fips-imm-init.csv',
  col_types = 'cnn'
)

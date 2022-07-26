cases  <- example_ct_data('cases')
hosp   <- example_ct_data('hosp')
RR     <- example_ct_data('RR')
boost  <- example_ct_data('boost')
deaths  <- example_ct_data('deaths')

imm_init <- get_imm_init("Connecticut")

config <- covidestim(nweeks = nrow(cases), region = 'Connecticut',
                  start_p_imm = imm_init$start_p_imm,
                  cum_p_inf_init = imm_init$cum_p_inf_init) +
  input_cases(cases) +
  input_hosp(hosp) + 
  input_rr(RR) + 
  input_boost(boost) +
  input_deaths(deaths)

result <- runOptimizer(config, cores=1, tries=1)
summarized <- summary(result)

test_that("Optimizer summary gets dates and data.available right", {

  expect_equal(1:4, which(!summarized$data_available))

  expect_equal(
    summarized$date[1],
    data.table::first(example_ct_data('cases')$date) -
      lubridate::days(7*4)
  )
})

test_that("Renaming of variables preserves value of variables", {
  expect_equal(
    summarized$deaths[1],
    result$result$deaths[1]
  )
})

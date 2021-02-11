config <- covidestim(
  ndays = 120,
  ndays_before = 28,
  region = 'New York',
  pop = get_pop('New York')
) + input_cases(example_nyc_data('cases')) +
    input_deaths(example_nyc_data('deaths'))

result <- runOptimizer(config, cores=1, tries=1)
summarized <- summary(result)

test_that("Optimizer summary gets dates and data.available right", {

  expect_equal(1:28, which(!summarized$data.available))

  expect_equal(
    summarized$date[1],
    first(example_nyc_data('cases')$date) -
      lubridate::days(28)
  )
})

test_that("Renaming of variables preserves value of variables", {
  expect_equal(
    summarized$deaths[1],
    result$result$new_die[1]
  )
})

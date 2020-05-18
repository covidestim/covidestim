test_that("bad parameters get thrown out", {

  outlandish_values <- list(0, 1, 2, -1, "hi", mtcars)

  purrr::walk(
    outlandish_values,
    ~expect_error(
      simulated_data(p_cases_die = .)
    )
  )

  purrr::walk(
    outlandish_values,
    ~expect_error(
      simulated_data(p_cases_diagnosed = .)
    )
  )

  purrr::walk(
    outlandish_values,
    ~expect_error(
      simulated_data(p_deaths_diagnosed = .)
    )
  )
})

test_that("there are no duplicate days in valid output", {
  d <- simulated_data()

  expect_equal(d$true$report_day,     unique(d$true$report_day))
  expect_equal(d$reported$report_day, unique(d$reported$report_day))
})

test_that("reporting data and truth data don't conflict", {
  d_original   <- simulated_data()
  d_cumulative <- d_original

  vs <- c("cases", "deaths")

  # Get cumsums for the days, but first, sort by reporting day
  d_cumulative$true     <- arrange(d_original$true,     report_day)
  d_cumulative$reported <- arrange(d_original$reported, report_day)
  d_cumulative$true     <- mutate_at(d_original$true,     vs, cumsum)
  d_cumulative$reported <- mutate_at(d_original$reported, vs, cumsum)

  # Join reported and truth data together
  joined <- dplyr::left_join(d_cumulative$true,
                             d_cumulative$reported,
                             by = "report_day",
                             suffix = c(".true", ".reported"))

  # Replace the NA values (introduced because the reporting lags the first
  # day of truth data) with 0s
  replace_na(
    joined,
    list(
      cases.true      = 0,
      deaths.true     = 0,
      cases.reported  = 0,
      deaths.reported = 0
    )
  ) -> joined

  # The cumulative number of reported deaths should always be l.t.e. the
  # number of true deaths, by that point. The same is true for cases
  expect_true(all(joined$deaths.reported <= joined$deaths.true))
  expect_true(all(joined$cases.reported  <= joined$cases.true))
})

test_that("reporting and truth data is returned already sorted", {
  d <- simulated_data()

  expect_identical(d$true,     arrange(d$true,     report_day))
  expect_identical(d$reported, arrange(d$reported, report_day))
})

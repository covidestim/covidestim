test_that("bad dates don't work", {
  library(lubridate)

  dates <- seq(now()-days(100), now(), by='1 day')[1:100]

  data.frame(
    date=format(dates, '%Y:%m:%d'),
    observations=rep(3, 100)
  ) -> d1

  data.frame(
    date=rep("10/10/2001", times=100),
    observations=rep(3,100)
  ) -> d2

  data.frame(
    date=rep("10/10/2001", times=100),
    observations=rep(3,100)
  ) -> d              

  expect_error(input_cases(d1))
  expect_error(input_cases(d2))
  expect_error(input_cases(d3))
})

test_that("duplicate dates don't work", {
  library(lubridate)

  dates <- seq(now()-days(100), now(), by='1 day')[1:100]
  dates <- c(dates[1:50], dates[50], dates[51:100])

  data.frame(
    date=dates,
    observations=rep(3, 101)
  ) -> d1

  expect_error(input_cases(d3))
})

test_that("example data validates", {

  d_cases  <- example_nyc_data("cases")
  d_deaths <- example_nyc_data("deaths")
  d_hosps  <- example_nyc_data("hosps")

  N_days <- nrow(d_cases)

  expect_silent(icas <- input_cases(d_cases))
  expect_silent(idth <- input_deaths(d_deaths))
  expect_silent(ihos <-input_hospitalizations(d_hosps))

  expect_silent(
    covidcast(N_days = N_days, N_days_delay = 10) +
      icas + idth + ihos
  )
})

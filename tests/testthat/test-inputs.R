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
  expect_error(input_cases(d1))
  expect_error(input_cases(d1))
  expect_error(input_deaths(d2))
  expect_error(input_deaths(d1))
  expect_error(input_deaths(d2))
  expect_error(input_fracpos(d3))
  expect_error(input_fracpos(d2))
  expect_error(input_fracpos(d3))
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
  expect_error(input_deaths(d3))
  expect_error(input_fracpos(d3))
})

test_that("example data validates", {

  d_cases  <- example_nyc_data("cases")
  d_deaths <- example_nyc_data("deaths")

  N_days <- nrow(d_cases)

  expect_silent(icas <- input_cases(d_cases))
  expect_silent(idth <- input_deaths(d_deaths))

  expect_silent(
    covidcast(N_days = N_days, N_days_before = 10) +
      icas + idth
  )
})

test_that("obs/rep causes change to underlying Stan configuration", {

  d_cases  <- example_nyc_data("cases")
  N_days <- nrow(d_cases)

  d_deaths <- example_nyc_data("deaths")
  d_fracpos <- tibble::tibble(
    date = seq(now() - N_days, now(), by = '1 day'),
    observation = runif(N_days, 0, 1)
  )


  # Configure things as reported
  expect_silent(icas     <- input_cases(d_cases))
  expect_silent(idth     <- input_deaths(d_deaths))
  expect_silent(ifracpos <- input_fracpos(d_fracpos))

  # Should show up as an attribute
  expect_equal(attr(icas, 'date_type'), "reported")
  expect_equal(attr(idth, 'date_type'), "reported")
  expect_equal(attr(ifracpos, 'date_type'), "reported")

  # Everything should succeed
  expect_silent(
    cfg <- covidcast(N_days = N_days, N_days_before = 10) + icas
  )

  # Everything should succeed
  expect_silent(
    cfg <- covidcast(N_days = N_days, N_days_before = 10) + idth
  )

  # Everything should succeed
  expect_silent(
    cfg <- covidcast(N_days = N_days, N_days_before = 10) + icas + idth
  )

  # Everything should succeed
  expect_silent(
    cfg <- covidcast(N_days = N_days, N_days_before = 10) + idth + icas
  )

  # You should be able to see the `*_rep` flags set now
  expect_equal(cfg$config$obs_cas_rep, 1)
  expect_equal(cfg$config$obs_die_rep, 1)

  # Everything should succeed
  expect_silent(
    cfg <- covidcast(N_days = N_days, N_days_before = 10) + ifracpos
  )

  expect_equal(cfg$config$frac_pos, c(rep(0, 10), d_fracpos$observation))

  #############################################################################

  # This next example checks to make sure that you can input "occurred" data
  # just as well as you can input "reported" data

  expect_silent(icas <- input_cases(d_cases,   type = "occurred"))
  expect_silent(idth <- input_deaths(d_deaths, type = "occurred"))

  # Should show up as an attribute
  expect_equal(attr(icas, 'date_type'), "occurred")
  expect_equal(attr(idth, 'date_type'), "occurred")

  # Everything should succeed
  expect_silent(
    cfg <- covidcast(N_days = N_days, N_days_before = 10) + icas + idth
  )

  # You should be able to see the `*_rep` flags UNset now
  expect_equal(cfg$config$obs_cas_rep, 0)
  expect_equal(cfg$config$obs_die_rep, 0)
})

test_that("bad `type` arguments don't validate to valid input objects", {

  d_cases  <- example_nyc_data("cases")
  d_deaths <- example_nyc_data("deaths")

  # Configure things as reported
  expect_error(icas <- input_cases(d_cases,  type = "LOL"))
  expect_error(idth <- input_deaths(d_deaths,type = "CAT"))
})

test_that("Not specifying fracpos results in the correct internal representation", {
  cfg <- covidcast(N_days = 50, N_days_before = 10)

  expect_equal(cfg$config$frac_pos, rep(0, 60))
  expect_null(cfg$config$frac_pos_user)
})

test_that("Weekend effect is represented correctly internally", {
  library(lubridate)
  library(magrittr)

  d_cases  <- example_nyc_data("cases")
  N_days <- nrow(d_cases)

  d_deaths <- example_nyc_data("deaths")

  # Configure things as reported
  expect_silent(icas     <- input_cases(d_cases))
  expect_silent(idth     <- input_deaths(d_deaths))

  # Everything should succeed
  expect_silent(
    cfg <- covidcast(N_days = N_days, N_days_before = 10) + icas + idth
  )

  first_day <- d_cases$date[1]
  first_week <- seq(first_day, first_day + days(6), by = '1 day')
  first_week_days <- wday(first_week)
  first_week_isweekend <- ifelse(first_week_days %in% c(6,7), 1, 0)

  expect_equal(first_week_isweekend, cfg$config$is_weekend[11:17])
})

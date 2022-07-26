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
})

test_that("example data validates", {

  d_cases  <- example_ct_data("cases")
  d_deaths <- example_ct_data("deaths")

  nweeks <- nrow(d_cases)
  imm_init <- get_imm_init("Connecticut")
  
  expect_silent(icas <- input_cases(d_cases))
  expect_silent(idth <- input_deaths(d_deaths))

  expect_silent(
    covidestim(nweeks = nrow(d_cases), region = 'Connecticut',
               start_p_imm = imm_init$start_p_imm,
               cum_p_inf_init = imm_init$cum_p_inf_init) +
      icas + idth
  )
})

test_that("obs/rep causes change to underlying Stan configuration", {

  d_cases  <- example_ct_data("cases")
  nweeks <- nrow(d_cases)
  imm_init <- get_imm_init("Connecticut")
  
  d_deaths <- example_ct_data("deaths")

  # Configure things as reported
  expect_silent(icas     <- input_cases(d_cases))
  expect_silent(idth     <- input_deaths(d_deaths))

  # Should show up as an attribute
  expect_equal(attr(icas, 'date_type'), "reported")
  expect_equal(attr(idth, 'date_type'), "reported")

  # Everything should succeed
  expect_silent(
    cfg <- covidestim(nweeks = nrow(d_cases), region = 'Connecticut',
               start_p_imm = imm_init$start_p_imm,
               cum_p_inf_init = imm_init$cum_p_inf_init) + icas
  )

  # Everything should succeed
  expect_silent(
    cfg <- covidestim(nweeks = nrow(d_cases), region = 'Connecticut',
                      start_p_imm = imm_init$start_p_imm,
                      cum_p_inf_init = imm_init$cum_p_inf_init) + idth
  )

  # Everything should succeed
  expect_silent(
    cfg <- covidestim(nweeks = nrow(d_cases), region = 'Connecticut',
                      start_p_imm = imm_init$start_p_imm,
                      cum_p_inf_init = imm_init$cum_p_inf_init) + icas + idth
  )

  # Everything should succeed
  expect_silent(
    cfg <- covidestim(nweeks = nrow(d_cases), region = 'Connecticut',
                      start_p_imm = imm_init$start_p_imm,
                      cum_p_inf_init = imm_init$cum_p_inf_init) + idth + icas
  )

  # You should be able to see the `*_rep` flags set now
  expect_equal(cfg$config$obs_cas_rep, 1)
  expect_equal(cfg$config$obs_die_rep, 1)

  # Everything should succeed
  expect_silent(
    cfg <- covidestim(nweeks = nrow(d_cases), region = 'Connecticut',
                      start_p_imm = imm_init$start_p_imm,
                      cum_p_inf_init = imm_init$cum_p_inf_init)
  )

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
    cfg <- covidestim(nweeks = nrow(d_cases), region = 'Connecticut',
                      start_p_imm = imm_init$start_p_imm,
                      cum_p_inf_init = imm_init$cum_p_inf_init) + icas + idth
  )

  # You should be able to see the `*_rep` flags UNset now
  expect_equal(cfg$config$obs_cas_rep, 0)
  expect_equal(cfg$config$obs_die_rep, 0)
})

test_that("bad `type` arguments don't validate to valid input objects", {

  d_cases  <- example_ct_data("cases")
  d_deaths <- example_ct_data("deaths")

  # Configure things as reported
  expect_error(icas <- input_cases(d_cases,  type = "LOL"))
  expect_error(idth <- input_deaths(d_deaths,type = "CAT"))
})


test_that("bad regions don't work", {
  expect_error(get_ifr_raw('Conn'))
  expect_error(get_ifr_raw('CT'))
  expect_error(get_ifr_raw('PR'))
  expect_error(get_ifr_raw('Puerto rico'))
  expect_error(get_ifr_raw('connecticut'))
  expect_error(get_ifr_raw('9009'))
  expect_error(get_ifr_raw(09009))
  expect_error(get_ifr_raw(09))
  expect_error(get_ifr_raw("09"))
})

test_that("good regions do work", {
  expect_length(get_ifr_raw('Connecticut'), 1)
  expect_length(get_ifr_raw('09009'), 1)
})

# test_that("start_date functions as expected", {
# 
#   startDate <- as.Date('2020-01-01')
# 
#   # Extent of the IFR data is 2020-01-01 - 2021-12-31
#   invalidDates <- as.Date(c('2019-12-31'), c('2022-01-01'))
# 
#   d <- get_ifr_raw('09009', startDate)
# 
#   expect_equal(min(d$date), startDate)
# 
#   expect_error(get_ifr_raw('09009', invalidDates[1]))
#   expect_error(get_ifr_raw('09009', invalidDates[2]))
# })

test_that('gen_ifr_adjustments outputs vectors of correct length', {

  N_weeks_before <- 5

  # 2020-01-22 is the earliest date of testing data for all counties in the
  # JHU dataset, and all states in the CTP dataset. Note that both of these
  # data have been cleaned, which *could*, but probably didn't, remove
  # days in the beginning.
  result <- gen_ifr_adjustments(as.Date('2020-01-22'), N_weeks_before, 'New York')

  expect_equal(length(result), 3)

  expect_equal(
    length(result$ifr_adj),
    seq.Date(
      as.Date('2020-01-22') - lubridate::weeks(N_weeks_before),
      as.Date('2022-12-31'),
      by = '1 week'
    ) %>% length
  )

  expect_equal(
    result$N_ifr_adj,
    seq.Date(
      as.Date('2020-01-22') - lubridate::weeks(N_weeks_before),
      as.Date('2022-12-31'),
      by = '1 week'
    ) %>% length
  )
})

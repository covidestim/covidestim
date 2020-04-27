library(lubridate)

test_that("bad dates don't work", {

  dates <- seq(now()-days(100), now(), by='1 day')

  d1 <- data.frame(date=format(dates, '%Y:%m:%d'), observs=rep(3, 101))

  expect_error(covidcast(N_days=101) + input_cases(d1))
})

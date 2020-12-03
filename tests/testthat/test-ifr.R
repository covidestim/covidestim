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
  expect_length(get_ifr_raw('Connecticut'), 2)
  expect_length(get_ifr_raw('09009'), 2)
})


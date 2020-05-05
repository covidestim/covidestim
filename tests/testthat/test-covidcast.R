test_that("addition is commutative", {
  d_deaths <- example_nyc_data('deaths')
  nd       <- nrow(d_deaths)

  expect_silent(
    form1 <- covidcast(N_days = nd) + 
      input_deaths(d_deaths) + 
      priors_progression(sym_prg_delay = c(5, 0.5))
  )

  expect_silent(
    form2 <- covidcast(N_days = nd) + 
      priors_progression(sym_prg_delay = c(5, 0.5)) +
      input_deaths(d_deaths)
  )

  expect_identical(form1, form2)
})

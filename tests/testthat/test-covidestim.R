test_that("addition is commutative", {
  d_deaths <- example_nyc_data('deaths')
  nd       <- nrow(d_deaths)

  expect_silent(
    form1 <- covidestim(ndays = nd) + 
      input_deaths(d_deaths) + 
      priors_progression(sym_prg_delay = c(5, 0.5))
  )

  expect_silent(
    form2 <- covidestim(ndays = nd) + 
      priors_progression(sym_prg_delay = c(5, 0.5)) +
      input_deaths(d_deaths)
  )

  # expect_mapequal(form1$config, form2$config)
})

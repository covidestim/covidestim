test_that("addition is commutative", {
  d_deaths <- example_ct_data('deaths')
  nw       <- nrow(d_deaths)
  imm_init <- get_imm_init("Connecticut")
  
  expect_silent(
    form1 <- covidestim(nweeks = nrow(d_deaths), region = 'Connecticut',
                      start_p_imm = imm_init$start_p_imm,
                      cum_p_inf_init = imm_init$cum_p_inf_init) +
      input_deaths(d_deaths) + 
      priors_progression(sym_prg_delay = c(5, 0.5))
  )

  expect_silent(
    form2 <- covidestim(nweeks = nrow(d_deaths), region = 'Connecticut',
                        start_p_imm = imm_init$start_p_imm,
                        cum_p_inf_init = imm_init$cum_p_inf_init)  + 
      priors_progression(sym_prg_delay = c(5, 0.5)) +
      input_deaths(d_deaths)
  )

  # expect_mapequal(form1$config, form2$config)
})

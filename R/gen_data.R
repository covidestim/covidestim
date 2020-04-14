set.seed(123)

genData <- function(diagData)
{
  datList <- list()

  list(
    N_conf_cases = length(diagData$diagnosis_day),

    #n days to model before first case
    N_days = max(diagData$diagnosis_day) + diagData$N_days_before,

    # model days after (for things like death, hospitalizations)
    N_days_extra = 1,

    # maximum value for delay distribution
    Max_delay = max(diagData$days_delay),

    cases_test_day = diagData$diagnosis_day + diagData$N_days_before,

    # NEED A PRIOR ON THIS
    cases_days_delay = diagData$days_delay,

    ## Priors Parameters of random walk in log space <- new infections per day
    # mean of log daily infections on day 1
    pri_log_new_inf_0_mu = -2,

    # sd of log daily infections on day 1
    pri_log_new_inf_0_sd = 2,

    # drift gives some direction to random walk. stronger prior here.
    # mean of daily change in log infections
    pri_log_new_inf_drift_mu = 0,
    # hyper prior on random walk
    # sd of daily change in log infections
    pri_log_new_inf_drift_sd = 1,
    # mean of daily change in log infections
    pri_sigma_deriv1_log_new_inf_sd = 0.5,

    # priors on the second derivative; penalizes sharp changes in random walk;
    # gives rw momentum
    pri_deriv2_log_new_inf_sd = 0.1,

    # probability of transitioning between states unlike ODE, not a constant
    # rate in transmission -- allows us to model delays probability of
    # transitioning into the next stage or recovery
    pri_p_sym_if_inf_mn = 0.69,
    # sum of the beta coeffecients
    pri_p_sym_if_inf_ss = 100,

    # upper bound from CDC paper
    pri_p_hos_if_sym_mn = 0.31,
    pri_p_hos_if_sym_ss = 100,

    # upper bound from CDC paper
    pri_p_die_if_hos_mn = 0.03,
    pri_p_die_if_hos_ss = 100,

    nb_yes = 0
  ) -> datList

  # Delay from infection to symptoms
  # https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported
  # Median incubation 5.1 days (4.5, 5.8 );
  # 97.5 percentile 11.5 days (8.2, 15.6);
  # 2.5 percentile 2.2 days (1.8, 2.9)

  # functions to get parameters from the Lauer paper
  exp(
    stats::optim(
      c(-1, -1),
      function(zz) {
        z   <- exp(zz)
        prd <- qgamma(c(1, 20, 39)/40, 1/z[1]/z[1], 1/z[1]/z[1]/z[2])
        sum( (prd - c(2.2, 5.1, 11.5))^2 * c(1, 5, 1) )
      },
      method="BFGS"
    )$par
  ) -> opt_par1

  # qgamma(c(1, 20, 39)/40, opt_par1[1]^-2, opt_par1[1]^-2/opt_par1[2])

  stats::optimize(
    function(z) {
      prd <- qgamma(c(1, 39)/40, 1/z[1]/z[1], 1/z[1]/z[1]/5.1)
      sum((prd - c(4.5, 5.8))^2)
    },
    c(0.01, 1)
  )$minimum -> opt_par2

  # qgamma(c(1, 39)/40, opt_par2^-2, opt_par2^-2/5.1)

  # gammapar <- function(tgt) {
  #   tgt <- as.numeric(tgt)
  #   mn <- tgt[1]
  #   cir <- (tgt[3] - tgt[2])
  #   xopt <- function(b, mn, cir) {
  #     cir2 <- qgamma(c(1, 39)/40, mn * b, b)
  #     cir2 <- cir2[2] - cir2[1]
  #     (cir2 - cir)^2
  #   }
  #   zz <- optimize(xopt, c(0.1, 1e+05), mn = mn, cir = cir)$minimum
  #   c(zz * mn, zz)
  # }

  # shape is fixed, mean can vary
  list(
    pri_inf_prg_delay_mn_mn = opt_par1[2],
    pri_inf_prg_delay_mn_cv = opt_par2,
    inf_prg_delay_cv        = opt_par1[1],

    # prior on the mean of the gamma distribution is distributed gamma with a
    # mean and a cv we could simplify by fixing the mean OR assuming time to
    # recovery and time to progression is same
    pri_sym_prg_delay_mn_mn = 4,
    pri_sym_prg_delay_mn_cv = 0.1,

    # cv on the gamma distribution
    sym_prg_delay_cv        = 0.5,

    pri_hos_prg_delay_mn_mn = 5,
    pri_hos_prg_delay_mn_cv = 0.1,
    hos_prg_delay_cv        = 0.5,

    pri_inf_res_delay_mn_mn = 14,
    pri_inf_res_delay_mn_cv = 0.1,
    inf_res_delay_cv        = 0.5,

    pri_sym_res_delay_mn_mn = 6,
    pri_sym_res_delay_mn_cv = 0.1,
    sym_res_delay_cv        = 0.5,

    pri_hos_res_delay_mn_mn = 6,
    pri_hos_res_delay_mn_cv = 0.1,
    hos_res_delay_cv        = 0.5,

    pri_report_delay_mn_mn  = 7,
    pri_report_delay_mn_cv  = 0.9,

    pri_report_delay_cv_mn  = 0.5,
    pri_report_delay_cv_cv  = 0.9,

    # daily probability of diagnosis as a function of individual health state
    # currently parameterizing beta distributions, but we need to add a
    # function of time probability of diagnosis is a logistic function with
    # intercepts by health state
    pri_p_diag_if_inf_mn    = 0.01,
    pri_p_diag_if_inf_ss    = 100,
    pri_p_diag_if_sym_mn    = 0.1,
    pri_p_diag_if_sym_ss    = 100,
    pri_p_diag_if_hos_mn    = 0.6,
    pri_p_diag_if_hos_ss    = 100
  ) -> moreParams

  datList <- purrr::splice(datList, moreParams)

  datList
}

defaultData <- function() genData(genFakeData())

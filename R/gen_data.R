set.seed(123)

#' @export
#' High level description of the function
#'
#' More extended description of the function
#'
#' @param param1 Description
#'
#' @param param2 Description
#'
#' @param param3 Description
#'
#' @return The return value
#'
#' @examples
#' print(mtcars)
#' @export
genData <- function(diagData)
{
  N_days <- max(diagData$diagnosis_day) + diagData$N_days_before
  N_days_extra <- 1

  # Number of parameters used to specify spline
  n_spl_par <- 10

  splines::bs(
    1:(N_days + N_days_extra),
    df        = n_spl_par,
    degree    = 3,
    intercept = T
  ) -> des_mat

  # this produces a cubic b-spline with n_spl_par basis functions
  spl_basis <- as.matrix(as.data.frame(des_mat))

  # The first set of components of 'datList'
  list(
    N_conf_cases = length(diagData$diagnosis_day),

    #n days to model before first case
    N_days = N_days,

    # model days after (for things like death, hospitalizations)
    N_days_extra = N_days_extra,

    # maximum value for delay distribution
    Max_delay = max(diagData$days_delay),

    cases_test_day = diagData$diagnosis_day + diagData$N_days_before,

    # NEED A PRIOR ON THIS
    cases_days_delay = diagData$days_delay,

    ## Priors Parameters of random walk in log space <- new infections per day
    # mean of log daily infections on day 1
    # pri_log_new_inf_0_mu = -2,
 
    # sd of log daily infections on day 1
    # pri_log_new_inf_0_sd = 2,

    # drift gives some direction to random walk. stronger prior here.
    # mean of daily change in log infections
    #pri_log_new_inf_drift_mu = 0,
    # hyper prior on random walk
    # sd of daily change in log infections
    # pri_log_new_inf_drift_sd = 1,
    # mean of daily change in log infections
    # pri_sigma_deriv1_log_new_inf_sd = 0.5,

    # priors on the second derivative; penalizes sharp changes in random walk;
    # gives rw momentum
    # pri_deriv2_log_new_inf_sd = 0.1,

    spl_basis = spl_basis,
    n_spl_par = n_spl_par,

    # transitions
    pri_p_sym_if_inf_a = 69,
    pri_p_sym_if_inf_b = 31,
    pri_p_hos_if_sym_a = 31,
    pri_p_hos_if_sym_b = 69,
    pri_p_die_if_hos_a = 3,
    pri_p_die_if_hos_b = 97,

    # poisson or negative binomial 
    nb_yes = 0,
    rw_yes = 0
  ) -> datList

  # The second set of components of 'datList'
  list(
    # delay to progression
    inf_prg_delay_shap = 5.202,
    inf_prg_delay_rate = 0.946,
    sym_prg_delay_shap = 51.47,
    sym_prg_delay_rate = 4.679,
    hos_prg_delay_shap = 91.64,
    hos_prg_delay_rate = 10.41,

    # delay to recovered  
      # inf to res not well supported in data, guess: 
      # assume mean 10 days, CI 5, 13
      # based on 5 days inf -> sym, 5 days sym -> res (CI 2,8)
      # adjust later? 
    inf_res_delay_shap = 23.83, 
    inf_res_delay_rate = 2.383,
    sym_res_delay_shap = 10.50,
    sym_res_delay_rate = 2.099,
    hos_res_delay_shap = 60.86,
    hos_res_delay_rate = 3.567,

    # report delay: assumed
    pri_report_delay_shap = 7,
    pri_report_delay_rate = 0.9,

    # probability of diagnosis 
      # assumed
    pri_p_diag_if_inf_a = 1,
    pri_p_diag_if_inf_b = 99,
    pri_p_diag_if_sym_a = 60,
    pri_p_diag_if_sym_b = 40,
    pri_p_diag_if_hos_a = 95,
    pri_p_diag_if_hos_b = 5

  ) -> moreParams

  datList <- purrr::splice(datList, moreParams)

  datList
}

#' High level description of the function
#'
#' More extended description of the function
#'
#' @param param1 Description
#'
#' @param param2 Description
#'
#' @param param3 Description
#'
#' @return The return value
#'
#' @examples
#' print(mtcars)
#' @export
defaultData <- function() genData(genFakeData())

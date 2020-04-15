# Alias for making assertions on code
att <- assertthat::assert_that

# Combines 'defaults' and 'args' together. If 'args'==list(), it is stripped
# from the result. Any duplicate keys are resolved by selecting the last key
# specified. The result has class 'modelconfig'
splice_class <- function(defaults, args, class)
  structure(
    rlang::dots_list(
      !!! defaults,
      !!! args,
      .homonyms     = 'last',
      .ignore_empty = 'trailing'),
    class=class)

# Example model configuration object. It can be modified through overloading
# the addition operator
modelconfig <- function(...) {
  config <- list(...)

  # Declare that this is a class
  structure(config, class='modelconfig')
}

# A list with arbitrary values, of class 'priors'
priors <- function(...) structure(list(...), class='priors')

# An overloaded addition operator, dispatched on the type of the lhs argument.
# Flips the order of the arguments and calls 'modelconfig_add' to enable type
# matching on the rhs operand.
"+.modelconfig" <- function(a, b) {
  # Dispatching on the type of 'b', which should be 'priors' for now.
  modelconfig_add(b, a)
}

modelconfig_add <- function(rightside, leftside) UseMethod('modelconfig_add') 

# Specialization for 'priors' classes. Splices in all the keys from a 'priors'
# obj into a 'modelconfig' (leftside) object
modelconfig_add.priors <- function(rightside, leftside)
  rlang::dots_list(!!!leftside, !!!rightside, .homonyms='last')
          
# BUG: the following description needs information about what the defaults
# were sourced from

#' Priors for transitions
#'
#' This function returns a keyed list of priors related to progression to
#' various health states. Called with no arguments, the default values are
#' returned. The following arguments can be passed to create different priors:
#'
#' \code{inf_prg_delay_shap}
#' \code{inf_prg_delay_rate}
#' \code{sym_prg_delay_shap}
#' \code{sym_prg_delay_rate}
#' \code{hos_prg_delay_shap}
#' \code{hos_prg_delay_rate}
#' 
#' Default values of these priors go here.
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' myData <- defaultData() + priors_progression(inf_prg_delay_shap = 3.5)
#' @export
priors_transitions <- function(...) {
  args <- list(...)

  list(
    pri_p_sym_if_inf_a = 69,
    pri_p_sym_if_inf_b = 31,
    pri_p_hos_if_sym_a = 31,
    pri_p_hos_if_sym_b = 69,
    pri_p_die_if_hos_a = 3,
    pri_p_die_if_hos_b = 97
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
}

# BUG: the following description needs information about what the defaults
# were sourced from

#' Priors on delay to progression
#'
#' This function returns a keyed list of priors related to progression to
#' the infectious state. Called with no arguments, the default values are
#' returned. The following arguments can be passed to create different priors:
#'
#' \code{inf_prg_delay_shap}
#' \code{inf_prg_delay_rate}
#' \code{sym_prg_delay_shap}
#' \code{sym_prg_delay_rate}
#' \code{hos_prg_delay_shap}
#' \code{hos_prg_delay_rate}
#' 
#' Default values of these priors 
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' myData <- defaultData() + priors_progression(inf_prg_delay_shap = 3.5)
#' @export
priors_progression <- function(...) {
  args <- list(...)

  list(
    inf_prg_delay_shap = 5.202,
    inf_prg_delay_rate = 0.946,
    sym_prg_delay_shap = 51.47,
    sym_prg_delay_rate = 4.679,
    hos_prg_delay_shap = 91.64,
    hos_prg_delay_rate = 10.41
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
}

# inf to res not well supported in data, guess:
# assume mean 10 days, CI 5, 13
# based on 5 days inf -> sym, 5 days sym -> res (CI 2,8)
# adjust later?
# BUG: ^ This needs to be integrated into the documentation!

#' Priors on delay to recovered
#'
#' This function returns a keyed list of priors related to recovery from
#' the infectious state. Called with no arguments, the default values are
#' returned. The following arguments can be passed to create different priors:
#'
#' \code{inf_res_delay_shap}
#' \code{inf_res_delay_rate}
#' \code{sym_res_delay_shap}
#' \code{sym_res_delay_rate}
#' \code{hos_res_delay_shap}
#' \code{hos_res_delay_rate}
#' 
#' Default values of these priors 
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' myData <- defaultData() + priors_recovery(inf_res_delay_shap = 20.12)
#' @export
priors_recovery <- function(...) {
  args <- list(...)

  list(
    inf_res_delay_shap = 23.83,
    inf_res_delay_rate = 2.383,
    sym_res_delay_shap = 10.50,
    sym_res_delay_rate = 2.099,
    hos_res_delay_shap = 60.86,
    hos_res_delay_rate = 3.567
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
}

#' Priors on reporting delay
#'
#' This function returns a keyed list of priors related to reporting delay.
#' Called with no arguments, the default values are returned. The following
#' arguments can be passed to create different priors:
#'
#' \code{pri_report_delay_shap}
#' \code{pri_report_delay_rate}
#'
#' Default values should be explained here, as well as constraints on custom
#' values. (Default values are assumed?)
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' myData <- defaultData() + priors_reporting_delay(pri_report_delay_shap = 6)
#' @export
priors_reporting_delay <- function(...) {
  args <- list(...)

  list(
    pri_report_delay_shap = 7,
    pri_report_delay_rate = 0.9
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
}

#' Priors on probability of diagnosis
#'
#' This function returns a keyed list of priors related to probability of
#' diagnosis.  Called with no arguments, the default values are returned. The
#' following arguments can be passed to create different priors:
#'
#' \code{pri_p_diag_if_inf_a}
#' \code{pri_p_diag_if_inf_b}
#' \code{pri_p_diag_if_sym_a}
#' \code{pri_p_diag_if_sym_b}
#' \code{pri_p_diag_if_hos_a}
#' \code{pri_p_diag_if_hos_b}
#'
#' Default values should be explained here, as well as constraints on custom
#' values. (Default values are assumed?)
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' myData <- defaultData() + priors_diagnosis(pri_p_diag_if_inf_a = 1.2)
#' @export
priors_diagnosis <- function(...) {
  args <- list(...)

  list(
    pri_p_diag_if_inf_a = 1,
    pri_p_diag_if_inf_b = 99,
    pri_p_diag_if_sym_a = 60,
    pri_p_diag_if_sym_b = 40,
    pri_p_diag_if_hos_a = 95,
    pri_p_diag_if_hos_b = 5
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
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
  config <- rlang::dots_list(
    .homonyms = "error", # Ensure that no keys are entered twice

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

    # poisson or negative binomial
    nb_yes = 0,
    rw_yes = 0,

    !!! priors_transitions(),
    !!! priors_progression(),
    !!! priors_recovery(),
    !!! priors_reporting_delay(),
    !!! priors_diagnosis()
  )

  structure(config, class='modelconfig')
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
defaultData <- genData(genFakeData())


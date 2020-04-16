#' Model configuration object. It can be modified through overloading
#' the addition operator
#' @export
modelconfig <- function(...) {
  config <- list(...)

  # Declare that this is a class
  structure(config, class='modelconfig')
}

modelconfig_add <- function(rightside, leftside) UseMethod('modelconfig_add') 

#' Specialization for 'priors' classes.
#'
#' Splices in all the keys from a 'priors' obj into a 'modelconfig' (leftside)
#' object.
#'
#' @export
modelconfig_add.priors <- function(rightside, leftside)
  rlang::dots_list(!!!leftside, !!!rightside, .homonyms='last')
          
#' An overloaded addition operator, dispatched on the type of the lhs argument.
#'
#' Flips the order of the arguments and calls 'modelconfig_add' to enable type
#' matching on the rhs operand.
#'
#' @export
"+.modelconfig" <- function(a, b) {
  # Dispatching on the type of 'b', which should be 'priors' for now.
  modelconfig_add(b, a)
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
defaultConfig <- function() genData(genFakeData())

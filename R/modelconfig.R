#' Model configuration object. It can be modified through overloading
#' the addition operator
modelconfig <- function(...) {
  config <- list(...)

  # Declare that this is a class
  structure(config, class='modelconfig')
}

#' @export
print.modelconfig <- function(cc) {
  print.priors(cc$config, .tab = TRUE)
  print.inputs(cc$config, .tab = TRUE)
}

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

modelconfig_add <- function(rightside, leftside) UseMethod('modelconfig_add') 

#' Specialization for 'priors' classes.
#'
#' Splices in all the keys from a 'priors' obj into a 'modelconfig' (leftside)
#' object.
modelconfig_add.priors <- function(rightside, leftside) {
  cfg <- rlang::dots_list(!!!leftside, !!!rightside, .homonyms='last')

  validate.modelconfig(cfg)

  structure(cfg, class = "modelconfig")
}

modelconfig_add.input <- function(rightside, leftside) {

  cfg  <- leftside
  d    <- rightside

  integer_keys <- c("obs_cas", "obs_die")
  keys         <- c(integer_keys, "frac_pos")

  # Create a list of any existing data
  preexisting_inputs_lst <- plyr::compact(cfg[integer_keys])
  preexisting_inputs     <- length(preexisting_inputs_lst) > 0

  att(length(d) == 1)
  att(names(d) %in% keys)

  att(!is.null(d[[1]]$observation))

  att(
    !names(d) %in% names(preexisting_inputs_lst),
    msg = glue("Input {names(d)} cannot be added twice to covidcast config")
  )

  att(
    length(unique(d$date)) == length(d$date),
    msg = "Only one observation per date may be entered for each input type"
  )
  
  att(
    length(d[[1]]$date) == cfg$N_days,
    msg = glue("Number of observations in {names(d)} ({length(d[[1]]$date)}) was not equal to N_days ({cfg$N_days})")
  )

  if (preexisting_inputs) {
    # At least one type of input has been added so far. Check that the start
    # dates are the same
    att(min(d[[1]]$date) == cfg$first_date)
  }

  # Update the first
  cfg$first_date  <- min(cfg$first_date, min(d[[1]]$date), na.rm=TRUE)

  # Update the is_weekend vector
  # local({
  #   first_date_Date <- as.Date(cfg$first_date, origin = '1970-01-01')
  #   seq(
  #     # First day in N_days_before
  #     first_date_Date - lubridate::days(cfg$N_days_before), 
  #     # Last day in N_days
  #     first_date_Date + lubridate::days(cfg$N_days - 1), 
  #     by = '1 day'
  #   ) -> entire_period

  #   days_of_week <- purrr::map_dbl(entire_period, lubridate::wday)

  #   # In lubridate, by default, 6 and 7 are Saturday and Sunday,
  #   # respectively
  #   ifelse(days_of_week %in% c(6,7), 1, 0)
  # }) -> cfg$is_weekend

  data_key <- names(d)
  data_type_key <- glue("{data_key}_rep")

  if (data_type_key %in% c("obs_cas_rep", "obs_die_rep")) {
    cfg[[data_key]] <- as.integer(d[[1]]$observation)

    cfg[[data_type_key]] <-
      switch(attr(d, "date_type"), "reported" = 1, "occurred" = 0)
  } else if (data_key == "frac_pos") {
    # frac_pos is handled differently because we store two copies: the "user"
    # copy which has length `N_days`, and the copy that will be passed to the
    # model, which has length `N_days_before + N_days`.
    cfg$frac_pos_user <- d[[1]]$observation
    cfg$frac_pos      <- c(rep(0, cfg$N_days_before), d[[1]]$observation)
  } else {
    stop(glue("{data_key} is not a valid input to a covidcast configuration"))
  }

  cfg

  structure(cfg, class = "modelconfig")
}

print.inputs <- function(cfg, .tab = FALSE) {

  t <- ifelse(.tab, '\t', '')

  frmtr <- function(d) format(length(d), width = 4, justify = 'centre')

  status_cases <-
    ifelse(is.null(cfg$obs_cas), '[ ❌ ]', glue('[{frmtr(cfg$obs_cas)}]'))
  status_deaths <-
    ifelse(is.null(cfg$obs_die), '[ ❌ ]', glue('[{frmtr(cfg$obs_die)}]'))
  status_fracpos <-
    ifelse(is.null(cfg$frac_pos_user), '[ ❌ ]', glue('[{frmtr(cfg$frac_pos_user)}]'))

'Inputs:

{t}{status_cases}\tCases
{t}{status_deaths}\tDeaths
{t}{status_fracpos}\tFraction positive
' -> msg

  cat(glue(msg))
}

validate.modelconfig <- function(cfg) {

  N_days <- cfg$N_days # For brevity

  case_reporting_mean  <- cfg$pri_cas_rep_delay_shap/cfg$pri_cas_rep_delay_rate
  death_reporting_mean <- cfg$pri_die_rep_delay_shap/cfg$pri_die_rep_delay_rate

'Mean case reporting delay (relative to total days of data) was longer than
expected.

Your case reporting delay (`cas_rep_delay`) has a mean of {case_reporting_mean}
days, whereas the total days of data, N_days, was {N_days}. Consider adjusting
or removing your custom prior.
' -> case_reporting_warning

'Mean case reporting delay (relative to total days of data) was longer than
expected.

Your case reporting delay (`die_rep_delay`) has a mean of {death_reporting_mean}
days, whereas the total days of data, `N_days`, was {N_days}. Consider adjusting
or removing your custom prior.
' -> death_reporting_warning

  # att_w(case_reporting_mean < N_days,  glue(case_reporting_warning))
  # att_w(death_reporting_mean < N_days, glue(death_reporting_warning))
}

  #############################################################################
## warning that was here is now obsolete. 
  
genData <- function(N_days, N_days_before = 21, rho = 1) #new default value
{
  att(rho > 0 && rho <= 1)

  # The first set of components of 'datList'
  config <- rlang::dots_list(
    .homonyms = "error", # Ensure that no keys are entered twice

    #n days of data to model 
    N_days = as.integer(N_days),

    #n day to model before start of data
    N_days_before = as.integer(N_days_before),
    
    #max delay to allow the model to consider. 50 is recommended. 
    Max_delay = 50, 

    rho = rho,

    is_weekend = rep(0, N_days_before + N_days),
    
    # vectors of event counts; default to 0 if no input
    obs_cas = NULL, # vector of int by date. should have 0s if no event that day
    obs_die = NULL, # vector of int by date. should have 0s if no event that day
    frac_pos = rep(0, N_days_before + N_days), # vector of int by date. default is 0
    frac_pos_user = NULL,

    # first day of data, as determined by looking at input data. This allows 
    # matching the above^ case data to specific dates.
    first_date = NA,

    ## Priors Parameters of random walk in log space <- new infections per day
    # mean of log daily infections on day 1
    pri_log_new_inf_0_mu = 0,
 
    # sd of log daily infections on day 1
    pri_log_new_inf_0_sd = 10,

    # drift gives some direction to random walk. stronger prior here.
    # mean of daily change in log infections
    #pri_log_new_inf_drift_mu = 0,
    # hyper prior on random walk
    # sd of daily change in log infections
    # pri_log_new_inf_drift_sd = 1,

    # priors on the second derivative; penalizes sharp changes in random walk;
    # gives rw momentum
    pri_deriv1_log_new_inf_sd = 0.5,
    pri_deriv2_log_new_inf_sd = 0.05,

    #spl_basis = spl_basis,
    #n_spl_par = n_spl_par,

    # poisson or negative binomial
    nb_yes      = as.integer(1),
    obs_cas_rep = as.integer(0),  # This ~means FALSE in stan
    obs_die_rep = as.integer(0),  # This ~means FALSE in stan
  )

  # Adding the priors in separately is neccessary in order to make checks
  # on the priors that depend on the value of 'config$N_days' run
  structure(config, class='modelconfig') +
    structure(
      rlang::dots_list(
        !!! priors_transitions(),
        !!! priors_progression(),
        !!! priors_reporting_delays(),
        !!! priors_diagnosis(),
        !!! priors_diagnosis_delays_scale()
      ),
      class = 'priors'
    )
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
defaultConfig <- function(...) genData(...)

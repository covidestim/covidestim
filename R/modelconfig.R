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

  d    <- rightside
  cfg  <- leftside
  keys <- c("obs_cas", "obs_die", "obs_hos")

  # Create a list of any existing data
  preexisting_inputs_lst <- plyr::compact(cfg[keys])
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

  # Need as.integer to make things play nicely with 'rstan'
  cfg[[names(d)]] <- as.integer(d[[1]]$observation)

  # Update the first
  cfg$first_date  <- min(cfg$first_date, min(d[[1]]$date), na.rm=TRUE)

  structure(cfg, class = "modelconfig")
}

print.inputs <- function(cfg, .tab = FALSE) {

  t <- ifelse(.tab, '\t', '')

  frmtr <- function(d) format(length(d), width = 4, justify = 'centre')

  status_cases <-
    ifelse(is.null(cfg$obs_cas), '[ ❌ ]', glue('[{frmtr(cfg$obs_cas)}]'))
  status_deaths <-
    ifelse(is.null(cfg$obs_die), '[ ❌ ]', glue('[{frmtr(cfg$obs_die)}]'))
  status_hospitalizations <-
    ifelse(is.null(cfg$obs_hos), '[ ❌ ]', glue('[{frmtr(cfg$obs_hos)}]'))

'Inputs:

{t}{status_cases} Cases
{t}{status_deaths} Deaths
{t}{status_hospitalizations} Hospitalizations

' -> msg

  cat(glue(msg))
}

validate.modelconfig <- function(cfg) {

  N_days <- cfg$N_days # For brevity

  case_reporting_mean  <- cfg$pri_cas_rep_delay_shap/cfg$pri_cas_rep_delay_rate
  hosp_reporting_mean  <- cfg$pri_hos_rep_delay_shap/cfg$pri_hos_rep_delay_rate
  death_reporting_mean <- cfg$pri_die_rep_delay_shap/cfg$pri_die_rep_delay_rate

'Mean case reporting delay (relative to total days of data) was longer than
expected.

Your case reporting delay (`cas_rep_delay`) has a mean of {case_reporting_mean}
days, whereas the total days of data, N_days, was {N_days}. Consider adjusting
or removing your custom prior.
' -> case_reporting_warning

'Mean hospitalization reporting delay (relative to total days of data) was
longer than expected.

Your hospitalization reporting delay (`hos_rep_delay`) has a mean of
{hosp_reporting_mean} days, whereas the total days of data, N_days, was
{N_days}. Consider adjusting or removing your custom prior.
' -> hosp_reporting_warning

'Mean case reporting delay (relative to total days of data) was longer than
expected.

Your case reporting delay (`die_rep_delay`) has a mean of {death_reporting_mean}
days, whereas the total days of data, `N_days`, was {N_days}. Consider adjusting
or removing your custom prior.
' -> death_reporting_warning

  att_w(case_reporting_mean < N_days,  glue(case_reporting_warning))
  att_w(hosp_reporting_mean < N_days,  glue(hosp_reporting_warning))
  att_w(death_reporting_mean < N_days, glue(death_reporting_warning))

  #############################################################################

'Mean probability of diagnosis was specified as being higher for asymptomatic
individuals (expressed as the prior `p_diag_if_inf`) than for symptomatic
individuals (expressed as the prior `p_diag_if_sym`).

mean(`p_diag_if_inf`) = {mean_diag_if_inf}
mean(`p_diag_if_sym`) = {mean_diag_if_sym}

Consider adjusting or removing your custom priors.
' -> diagnosis_symptomatic_warning

'Mean probability of diagnosis was specified as being higher for symptomatic,
non-hospitalized individuals (expressed as the prior `p_diag_if_sym`) than for
hospitalized individuals (expressed as the prior `p_diag_if_hos`).

mean(`p_diag_if_sym`) = {mean_diag_if_sym}
mean(`p_diag_if_hos`) = {mean_diag_if_hos}

Consider adjusting or removing your custom priors.
' -> diagnosis_hospitalized_warning

  mean_diag_if_inf <- cfg$pri_p_diag_if_inf_a  / (cfg$pri_p_diag_if_inf_a  + cfg$pri_p_diag_if_inf_b)
  mean_diag_if_sym <- cfg$pri_p_diag_if_sym_a / (cfg$pri_p_diag_if_sym_a + cfg$pri_p_diag_if_sym_b)
  mean_diag_if_hos <- cfg$pri_p_diag_if_hos_a / (cfg$pri_p_diag_if_hos_a + cfg$pri_p_diag_if_hos_b)

  att_w(mean_diag_if_inf <= mean_diag_if_sym, glue(diagnosis_symptomatic_warning))
  att_w(mean_diag_if_sym <= mean_diag_if_hos, glue(diagnosis_hospitalized_warning))
}

genData <- function(N_days, N_days_delay = 10)
{
  # The first set of components of 'datList'
  config <- rlang::dots_list(
    .homonyms = "error", # Ensure that no keys are entered twice

    #n days of data to model 
    N_days = as.integer(N_days),

    #n day to model before start of data
    N_days_delay = as.integer(N_days_delay),
    
    # vectors of event counts; default to 0 if no input
    obs_cas = NULL, # vector of int by date. should have 0s if no event that day
    obs_hos = NULL, # vector of int by date. should have 0s if no event that day
    obs_die = NULL, # vector of int by date. should have 0s if no event that day

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
    obs_cas_rep = as.integer(0), 
    obs_hos_rep = as.integer(0), 
    obs_die_rep = as.integer(0),
  )

  # Adding the priors in separately is neccessary in order to make checks
  # on the priors that depend on the value of 'config$N_days' run
  structure(config, class='modelconfig') +
    structure(
      rlang::dots_list(
        !!! priors_transitions(),
        !!! priors_progression(),
        !!! priors_recovery(),
        !!! priors_reporting_delay(),
        !!! priors_diagnosis(),
        !!! priors_reporting_delays_new()
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

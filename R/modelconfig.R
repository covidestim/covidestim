#' Model configuration object. It can be modified through overloading
#' the addition operator
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
modelconfig_add.priors <- function(rightside, leftside) {
  cfg <- rlang::dots_list(!!!leftside, !!!rightside, .homonyms='last')

  validate.modelconfig(cfg)

  cfg
}

print.inputs <- function(cfg, .tab = FALSE) {

  t <- ifelse(.tab, '\t', '')

  status_cases <-
    ifelse(is.null(cfg$obs_cas), '[      ]', glue('[loaded]\t[{length(cfg$obs_cas)} obs]'))
  status_deaths <-
    ifelse(is.null(cfg$obs_die), '[      ]', glue('[loaded]\t[{length(cfg$obs_cas)} obs]'))
  status_hospitalizations <-
    ifelse(is.null(cfg$obs_hos), '[      ]', glue('[loaded]\t[{length(cfg$obs_cas)} obs]'))

'Inputs:

{t}Cases\t{status_cases}
{t}Deaths\t{status_deaths}
{t}Hospitalizations\t{status_hospitalizations}

' -> msg

  cat(glue(msg))
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

  cfg[[names(d)]] <- d[[1]]$observation
  cfg$first_date  <- min(cfg$first_date, min(d[[1]]$date), na.rm=TRUE)

  cfg
}

validate.modelconfig <- function(cfg) {

  N_days <- cfg$N_days # For brevity

  att_w(cfg$pri_inf_prg_delay_shap/cfg$pri_inf_prg_delay_rate < N_days/2,
   "Mean delay from infection to symptom onset (relative to total days of data) is longer than expected.")
  att_w(cfg$pri_sym_prg_delay_shap/cfg$pri_sym_prg_delay_rate < N_days/2,
   "Mean delay from symptom onset to hospitalization (relative to total days) of data is longer than expected.")
  att_w(cfg$pri_hos_prg_delay_shap/cfg$pri_hos_prg_delay_rate < N_days/1.5,
   "Mean delay from hospitalization to death (relative to total days of data) is longer than expected.")
  att_w(cfg$pri_inf_res_delay_shap/cfg$pri_inf_res_delay_rate < N_days/2,
   "Mean delay from asymptomatic infection to recovery (relative to total days of data) is longer than expected.")
  att_w(cfg$pri_sym_res_delay_shap/cfg$pri_sym_res_delay_rate < N_days/2,
   "Mean delay from symptom onset to recovery (relative to total days of data) is longer than expected.")
  att_w(cfg$pri_hos_res_delay_shap/cfg$pri_hos_res_delay_rate < N_days/1.5,
   "Mean delay from hospitalization to recovery (relative to total days of data) is longer than expected.")
  att_w(cfg$pri_cas_rep_delay_shap/cfg$pri_cas_rep_delay_rate < N_days/3,
   "Mean reporting delay for cases (relative to total days of data) is longer than expected.")
  att_w(cfg$pri_hos_rep_delay_shap/cfg$pri_hos_rep_delay_rate < N_days/3,
   "Mean reporting delay for hospitalizations (relative to total days of data) is longer than expected.")
  att_w(cfg$pri_die_rep_delay_shap/cfg$pri_die_rep_delay_rate < N_days/3,
   "Mean reporting delay for deaths (relative to total days of data) is longer than expected.")

  mean_sym_if_inf  <- cfg$pri_p_sym_if_inf_a  / (cfg$pri_p_sym_if_inf_a  + cfg$pri_p_sym_if_inf_b)
  mean_diag_if_sym <- cfg$pri_p_diag_if_sym_a / (cfg$pri_p_diag_if_sym_a + cfg$pri_p_diag_if_sym_b)
  mean_diag_if_hos <- cfg$pri_p_diag_if_hos_a / (cfg$pri_p_diag_if_hos_a + cfg$pri_p_diag_if_hos_b)

  att_w(mean_sym_if_inf <= mean_diag_if_sym,
        "Mean probability of diagnosis is higher for asymptomatic individuals than for symptomatic individuals.")
  att_w(mean_diag_if_sym <= mean_diag_if_hos,
        "Mean probability of diagnosis is higher for non-hospitalized individuals than for hospitalized individuals.")
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

genData <- function(diagData)
{
  N_days <- max(diagData$diagnosis_day) + diagData$N_days_before
  N_days_extra <- 1

  # The first set of components of 'datList'
  config <- rlang::dots_list(
    .homonyms = "error", # Ensure that no keys are entered twice

    #n days of data to model 
    N_days = N_days,

    #n day to model before start of data
    N_days_delay = 10,
    
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
    # mean of daily change in log infections
    pri_sigma_deriv1_log_new_inf_sd = 1,

    # priors on the second derivative; penalizes sharp changes in random walk;
    # gives rw momentum
    pri_deriv1_log_new_inf_sd = 0.5,
    pri_deriv2_log_new_inf_sd = 0.05,

    #spl_basis = spl_basis,
    #n_spl_par = n_spl_par,

    # poisson or negative binomial
    nb_yes = 0,
    obs_cas_rep = 0, 
    obs_hos_rep = 0, 
    obs_die_rep = 0,
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
        !!! priors_fixed()
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
defaultConfig <- function() genData(genFakeData())

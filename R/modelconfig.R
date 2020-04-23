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
modelconfig_add.priors <- function(rightside, leftside) {
  cfg <- rlang::dots_list(!!!leftside, !!!rightside, .homonyms='last')

  validate.modelconfig(cfg)

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

  # The first set of components of 'datList'
  config <- rlang::dots_list(
    .homonyms = "error", # Ensure that no keys are entered twice

    #n days of data to model 
    N_days = N_days,

    #n day to model before start of data
    N_days_delay = 10,
    
    # vectors of event counts; default to 0 if no input
    obs_cas = rep(0, N_days), # vector of int by date. should have 0s if no event that day
    obs_hos = rep(0, N_days), # vector of int by date. should have 0s if no event that day
    obs_die = rep(0, N_days), # vector of int by date. should have 0s if no event that day

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
#' @export
defaultConfig <- function() genData(genFakeData())

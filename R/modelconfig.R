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
    msg = glue("Input {names(d)} cannot be added twice to covidestim config")
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

  # Update the is_weekend vector, IFF `weekend = TRUE` was passed
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

  #   # In lubridate, by default, 7 and 1 are Saturday and Sunday,
  #   # respectively
  #   ifelse(days_of_week %in% c(7,1), 1, 0)
  # }) -> weekend_modifiers
  # 
  # if (cfg$user_weekend_effect)
  #   cfg$is_weekend <- weekend_modifiers

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
    stop(glue("{data_key} is not a valid input to a covidestim configuration"))
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

genData <- function(N_days, N_days_before = 28) #new default value
{

   N_days_ <- as.integer(N_days)
   N_days_before_ <- as.integer(N_days_before) - 7

  # The first set of components of 'datList'
  config <- rlang::dots_list(
    .homonyms = "error", # Ensure that no keys are entered twice
    
    # Data
    # vector of non-negative integer daily diagnoses, with a week of leading zeros
    obs_cas = NULL,   

    # vector of non-negative integer daily deaths, with a week of leading zeros
    obs_die = NULL,   

    # non-negative integer of data length  
    N_days = N_days_,

    # non-negative integer, length of burn-in, minus the week of leading zeros
    N_days_before = N_days_before_,

    # non-negative integer
    Max_delay = 6*7,      

  # Incidence time series
    #  real
    pri_log_new_inf_0_mu = 0, 

    # non-negative real
    pri_log_new_inf_0_sd = 10,

    # priors_list$serial_intvl_lnorm_pars[1], 
    pri_serial_i_a = 1.753997, 

    # priors_list$serial_intvl_lnorm_pars[2],
    pri_serial_i_b = 0.08787308, 

    # intercept for logRt
    pri_logRt_mu = 0, 

    # intercept for logRt
    pri_logRt_sd = 1.5, 

    # penalizes changes in Rt level
    pri_deriv1_spl_par_sd = 0.5, 

    # penalizes changes in Rt curvature
    pri_deriv2_spl_par_sd = 0.1, 

    # imported infections 
    pri_inf_imported_mu = 0,   

    # imported infections mean 0.5 per day
    pri_inf_imported_sd = 0.5/0.798,   

  # Transition probabilities  
    # priors_list$pct_inf_sym_beta_prior[1],  # positive real
    pri_p_sym_if_inf_a = 5.14303, 

    # priors_list$pct_inf_sym_beta_prior[2],  # positive real
    pri_p_sym_if_inf_b = 3.535966, 

    # priors_list$pct_sym_sev_beta_prior[1] , # positive real
    pri_p_sev_if_sym_a = 1.885352, 

    # priors_list$pct_sym_sev_beta_prior[2] , # positive real
    pri_p_sev_if_sym_b = 20.00246, 

    # priors_list$pct_sev_die_beta_prior[1] , # positive real
    pri_p_die_if_sev_a = 28.23883, 

    # priors_list$pct_sev_die_beta_prior[2] , # positive real
    pri_p_die_if_sev_b = 162.303, 

    # priors_list$ifr_beta_prior[1],          # positive real
    pri_p_die_if_inf_a = 72.14682, 

    # priors_list$ifr_beta_prior[2],          # positive real
    pri_p_die_if_inf_b = 10976.78, 

  # Probability of diagnosis (terms for beta priors)
    # positive real
    pri_p_diag_if_sev_a = 5,          

    # positive real
    pri_p_diag_if_sev_b = 2,          

    # positive real
    pri_rr_diag_sym_vs_sev_a = 2,     

    # positive real
    pri_rr_diag_sym_vs_sev_b = 2,      

  # Progression delays
    # priors_list$inf_prg_delay_gamma_pars[1], #  positive real
    inf_prg_delay_shap = 3.4131301, 

    # priors_list$inf_prg_delay_gamma_pars[2], # positive real
    inf_prg_delay_rate = 0.6051395, 

    # priors_list$sym_prg_delay_gamma_pars[1], # positive real
    sym_prg_delay_shap = 1.6236462, 

    # priors_list$sym_prg_delay_gamma_pars[1], # positive real
    sym_prg_delay_rate = 0.2175309, 

    # priors_list$sev_prg_delay_gamma_pars[1], # positive real
    sev_prg_delay_shap = 2.061454,  

    # priors_list$sev_prg_delay_gamma_pars[2], # positive real
    sev_prg_delay_rate = 0.227708,  

  # Scale diagnosis delay to progression delay (terms for beta priors)
    scale_dx_delay_sym_a = 2,  #  positive real
    scale_dx_delay_sym_b = 2,  #  positive real
    scale_dx_delay_sev_a = 2,  #  positive real
    scale_dx_delay_sev_b = 2,  #  positive real
  # Reporting delays (terms for lognormal distribution of beta prior terms)  
    pri_cas_rep_delay_shap_a = log(2.2), #  positive real
    pri_cas_rep_delay_rate_a = log(1),   #  positive real 
    pri_die_rep_delay_shap_a = log(2.2), #  positive real
    pri_die_rep_delay_rate_a = log(1),   #  positive real
    pri_cas_rep_delay_shap_b = 0.5,      #  positive real
    pri_cas_rep_delay_rate_b = 0.5,      #  positive real
    pri_die_rep_delay_shap_b = 0.5,      #  positive real
    pri_die_rep_delay_rate_b = 0.5,      #  positive real
  # Prior on phi  
    pri_inv_sqrt_phi = 1.0, #  positive real
  # Which data included  
    die_yes = 1,  #  0 or 1
    cas_yes = 1,  #  0 or 1
  # Are data by date of case report (1) or by date of test (0)?
    obs_cas_rep = 1, #  0 or 1
    obs_die_rep = 1, #  0 or 1
  # No. days of moving average 
    n_day_av = 5,   #  integer in 1:7 )

    !!! local({
      # Spline parameters for rt
      days_per_spl_par <- 4  # non-negative integer
      tot_days <- N_days_ + N_days_before_
      no_spl_par_rt <- ceiling(tot_days/days_per_spl_par)

      splines::bs(
        1:(days_per_spl_par*no_spl_par_rt),
        df = no_spl_par_rt,
        degree = 3,
        intercept = T
      ) -> des_mat

      list(
        N_spl_par = no_spl_par_rt,
        spl_basis = as.matrix(as.data.frame(des_mat))[1:tot_days,]
      )
    }),

    !!! local({
      tot_days <- N_days_ + N_days_before_
      N_spl_par_dx <- 6

      splines::bs(1:tot_days,
         df = N_spl_par_dx,
         degree = 3,
         intercept = T
      ) -> des_mat_dx  

      list(
        spl_basis_dx = as.matrix(as.data.frame(des_mat_dx)),
        N_spl_par_dx = N_spl_par_dx
      )
    })
  )

  structure(config, class='modelconfig')

  # Adding the priors in separately is neccessary in order to make checks
  # on the priors that depend on the value of 'config$N_days' run
  # structure(config, class='modelconfig') +
  #   structure(
  #     rlang::dots_list(
  #       !!! priors_transitions(),
  #       !!! priors_progression(),
  #       !!! priors_reporting_delays(),
  #       !!! priors_diagnosis(),
  #       !!! priors_diagnosis_delays_scale()
  #     ),
  #     class = 'priors'
  #   )
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

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

  integer_keys <- c("obs_cas", "obs_die", 
                    "obs_boost", "obs_hosp",
                    "ifr_vac_adj")
  keys         <- integer_keys

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
    # length(d[[1]]$date) == cfg$N_days,
    # msg = glue("Number of observations in {names(d)} ({length(d[[1]]$date)}) was not equal to N_days ({cfg$N_days})")
    length(d[[1]]$date) == cfg$N_weeks,
    msg = glue("Number of observations in {names(d)} ({length(d[[1]]$date)}) was not equal to N_weeks ({cfg$N_weeks})")
  )

  if (preexisting_inputs) {
    # At least one type of input has been added so far. Check that the start
    # dates are the same
    att(min(d[[1]]$date) == cfg$first_date)
  }

  # Update the first
  cfg$first_date  <- min(cfg$first_date, min(d[[1]]$date), na.rm=TRUE)
  
  if("lastDeathDate" %in% names(attributes(d))){
    lastDate <- attr(d, "lastDeathDate")
    diff <- lastDate - as.Date(cfg$first_date, origin = '1970-01-01')
    cfg$lastDeathWeek <- as.numeric(diff)%/%7
  } 

    if("lastHospDate" %in% names(attributes(d))){
    lastDate <- attr(d, "lastHospDate")
    diff <- lastDate - as.Date(cfg$first_date, origin = '1970-01-01')
    cfg$lastHospWeek <- as.numeric(diff)%/%7
  } 

  if("lastCaseDate" %in% names(attributes(d))){
    lastDate <- attr(d, "lastCaseDate")
    diff <- lastDate - as.Date(cfg$first_date, origin = '1970-01-01')
    cfg$lastCaseWeek <- as.numeric(diff)%/%7
  } 

  data_key <- names(d)
  data_type_key <- glue("{data_key}_rep")

  if (data_type_key %in% c("obs_cas_rep", "obs_die_rep", 
                           "obs_boost_rep", "obs_hosp_rep", 
                           "ifr_vac_adj_rep")) {
    if(data_type_key == "ifr_vac_adj_rep"){
      cfg[[data_key]] <- d[[1]]$observation 
    } else {
      cfg[[data_key]] <- as.integer(d[[1]]$observation)
    }

    cfg[[data_type_key]] <-
      switch(attr(d, "date_type"), "reported" = 1, "occurred" = 0)
  } else {
    stop(glue("{data_key} is not a valid input to a covidestim configuration"))
  }

  # The next several lines add state/county-specific IFR data to the
  # configuration passed to Stan. These lines have to be defined inside 
  # this function because they need to know the start date of the input data,
  # which is not known until the user adds `input_cases()/input_deaths()` to
  # the model configuration.
  ifr_adjustments <- gen_ifr_adjustments(
    first_date    = cfg$first_date %>% as.Date(origin = '1970-01-01'),
    N_weeks_before = cfg$N_weeks_before,
    region        = cfg$region
    )

  # Assign the results of the call to the `cfg` object
  cfg$ifr_adj       = ifr_adjustments$ifr_adj
  cfg$ifr_omi       = ifr_adjustments$ifr_omi
  cfg$ifr_adj_fixed = ifr_adjustments$ifr_adj_fixed
  cfg$N_ifr_adj     = ifr_adjustments$N_ifr_adj
  
  # pre-fill the ifr_vaccine_adjustment with ones, such that
  # the length of ifr_vac_adj matches N_ifr_adj
  if (data_key == "ifr_vac_adj"){
  # prefill <- rep(1, cfg$N_days_before)
  prefill <- rep(1, cfg$N_weeks_before)
  cfg$ifr_vac_adj   = c(prefill, cfg$ifr_vac_adj)
  }

  structure(cfg, class = "modelconfig")
}

print.inputs <- function(cfg, .tab = FALSE) {

  t <- ifelse(.tab, '\t', '')

  frmtr <- function(d) format(length(d), width = 4, justify = 'centre')

  status_cases <-
    ifelse(is.null(cfg$obs_cas), '[ x ]', glue('[{frmtr(cfg$obs_cas)}]'))
  status_deaths <-
    ifelse(is.null(cfg$obs_die), '[ x ]', glue('[{frmtr(cfg$obs_die)}]'))
  status_boost <-
    ifelse(is.null(cfg$obs_boost), '[ x ]', glue('[{frmtr(cfg$obs_boost)}]'))
  status_hospi <-
    ifelse(is.null(cfg$obs_hosp), '[ x ]', glue('[{frmtr(cfg$obs_hosp)}]'))
  status_rr <-
    ifelse(is.null(cfg$ifr_vac_adj), '[ x ]', glue('[{frmtr(cfg$ifr_vac_adj)}]'))

'Inputs:

{t}{status_cases}\tCases
{t}{status_deaths}\tDeaths
{t}{status_boost}\tBooster
{t}{status_hospi}\tHospitalizations
{t}{status_rr}\tRRvaccines
' -> msg

  cat(glue(msg))
}

validate.modelconfig <- function(cfg) {

  if (!is.character(cfg$region) || length(cfg$region) != 1)
    stop(glue::glue('You passed a `region` that is not a string!
                    Regions must be FIPS codes or state names'))

}

genData <- function(N_weeks, N_weeks_before = 28/7,
                    pop_size = 1e12, #new default value
                    n_spl_rt_knotwidth = 2, 
                    region,
                    cum_p_inf_init,
                    start_p_imm,
                    input_df
                    )
{
  n_spl_par_rt <- max(4,ceiling((N_weeks + N_weeks_before)/n_spl_rt_knotwidth))
  des_mat_rt <- splines::bs(1:(N_weeks + N_weeks_before), 
                            df=n_spl_par_rt, degree=3, intercept=T)

  n_spl_par_dx <- max(4,ceiling((N_weeks + N_weeks_before)/3)) 
  des_mat_dx <- splines::bs(1:(N_weeks + N_weeks_before), 
                         df=n_spl_par_dx, degree=3, intercept=T) 
  
  # The first set of components of 'datList'
  config <- rlang::dots_list(
    .homonyms = "error", # Ensure that no keys are entered twice

    #n days of data to model 
    # N_days = as.integer(N_days),
    N_weeks = as.integer(N_weeks),
    
    # Region being modeled. Must be a state name, or a FIPS code, passed as
    # a string.
    region = region,
    start_p_imm = start_p_imm,
    cum_p_inf_init = cum_p_inf_init,

    # These next three variables are NULLed because they aren't defined here.
    # They get defined when an input_*() is added, because it's there that
    # we first see what the user's first day of data is.
    ifr_adj       = NULL,
    ifr_omi       = NULL,
    ifr_adj_fixed = NULL,
    N_ifr_adj     = NULL,

    #n days to model before start of data
    # N_days_before = as.integer(N_days_before),
    N_weeks_before = as.integer(N_weeks_before),
    #max delay to allow the model to consider. 30 days is recommended; 5 weeks 
    Max_delay = ceiling(30/7), 

    # moving average for likelihood function 
    # N_days_av = N_days_av,
    
    # Whether to assume no reported cases and deats before the data (0 = no, 1 = yes) 
    pre_period_zero = 1,

    # Add population size to constrain susceptible population, large default assumes no constraint
    pop_size = pop_size, 
    
    # vectors of event counts; default to 0 if no input
    obs_cas = NULL, # vector of int by date. should have 0s if no event that day
    obs_hosp = NULL, # vector of int by date. should have 0s if no event that day
    obs_die = NULL, # vector of int by date. should have 0s if no event that day
    obs_boost = NULL, # vector of int by date. should have 0s if no event that day
    # the ifr_vaccine adjustment data
    ifr_vac_adj = NULL,
    # first day of data, as determined by looking at input data. This allows 
    # matching the above^ case data to specific dates.
    first_date = NA,
    # last death date and hosp date are initiated here; gets updated as deaths data gets added
    lastDeathWeek = N_weeks,
    lastHospWeek  = N_weeks,
    lastCaseWeek  = N_weeks,

    
    # Rt and new infections
    pri_log_infections_0_mu = 0,
    pri_log_infections_0_sd = 10,
    pri_logRt_mu = 0, # intercept for logRt
    pri_logRt_sd = 3.0, # sd of intercept for logRt
    pri_inf_imported_mu = 0,   # imported infections
    pri_inf_imported_sd = 0.5/0.798,   # imported infections 
    pri_deriv1_spl_par_sd = 0.5, # penalizes changes in Rt level
    pri_deriv2_spl_par_sd = 0.1, # penalizes changes in Rt curvature
    
    N_spl_par_rt = n_spl_par_rt, 
    spl_basis_rt = as.matrix(as.data.frame(des_mat_rt)),
    N_spl_par_dx = n_spl_par_dx, 
    spl_basis_dx = as.matrix(as.data.frame(des_mat_dx)),

    # indicates whether case or death data are being used 
    cas_yes = as.integer(1), 
    hosp_yes = as.integer(1), 
    die_yes = as.integer(1), 
    boost_yes = as.integer(1), 
    
    obs_cas_rep = as.integer(0),      # This ~means FALSE in stan
    obs_hosp_rep = as.integer(0),     # This ~means FALSE in stan
    obs_die_rep = as.integer(0),      # This ~means FALSE in stan
    obs_boost_rep = as.integer(0),      # This ~means FALSE in stan
    ifr_vac_adj_rep = as.integer(0),  # This ~means FALSE in stan
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

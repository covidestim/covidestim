# Equations for converting a mean and variance to the shape/rate parameters
# of a gamma distribution.

#' Convert mean/variance stats to shape/rate parameters of gamma distribution
#'
#' Using a method of moments approach, converts a sample with a specified
#' mean and variance to a parameterization of a gamma distribution
#'
#' @param mv A two-element numeric vector, where \code{[1]} represents mean
#'   and \code{[2]} represents variance.
#'
#' @return A number representing the shape or rate parameter
#' @export
pri_gamma_shape <- function(mv) (mv[1]^2)/mv[2]

#' @rdname pri_gamma_shape
pri_gamma_rate  <- function(mv) mv[1]/mv[2]

# A list with arbitrary values, of class 'priors'
priors <- function(...) structure(list(...), class='priors')

print.priors <- function(ps, .tab = FALSE) {
'Priors:

' -> msg

  cat(msg)

  for (idx in names(ps)) {
    if (stringr::str_detect(idx, '^pri') ||
        # Deal with "priors" that don't begin with "pri_"
        any(stringr::str_detect(idx, c("_prg_", "_res_")))) {

      # Random stuff for extra pretty-printing
      idx_better <- stringr::str_replace(idx, '^pri_', '')
      idx_better <- stringr::str_replace(idx_better, '_a', '\t[alpha]')
      idx_better <- stringr::str_replace(idx_better, '_b', '\t[beta]')
      idx_better <- stringr::str_replace(idx_better, '_shap', '\t[shape]')
      idx_better <- stringr::str_replace(idx_better, '_rate', '\t[rate]')

      if (.tab)
        cat('\t')

      cat(glue("{idx_better}\t{ps[[idx]]}\n\n"))
    }
  }
  cat('\n')
}

# Custom test. May experiment with 'on_failure(f) <- ' later on.
is_nonNegativeReal <- function(x, .element_names = NULL) {
  att(is.numeric(x))
  
  all(x >= 0)
}

# assertthat::on_failure(is_nonNegativeReal) <- function(call, env) {
#   idxs <- which(deparse(call$x) < 0)
#   
#   if (!is.null(env$.element_names))
#     els <- env$.element_names[idxs]
# 
# #  glue("{els} were not non-negative real numbers: {idxs}")
#   names(env)
# }

#' Helper function for manipulation of prior-related list-structures
#'
#' Assumes each argument is a two-element vector. Captures the name of each
#' argument. Creates two list-elements for each argument, using the 
#' naming convention given by '.postfix'. Prefixes these list-element keys
#' with the string '.prefix'. Returns the new list.
#'
#' @importFrom magrittr %>%
build_priors <- function(..., .postfix = c("_a", "_b"), .prefix = "") {

  att(is.character(.postfix))
  att(assertthat::is.string(.prefix))

  arg_names <- purrr::map_chr(rlang::enquos(...), rlang::quo_name)

  args_rekeyed <- stats::setNames(list(...), arg_names)

  purrr::lmap(
    args_rekeyed,
    ~rlang::dots_list(
      !! glue("{.prefix}{names(.)}{.postfix[1]}") := .[[1]][1],
      !! glue("{.prefix}{names(.)}{.postfix[2]}") := .[[1]][2],
      .homonyms = 'error'
    )
  )
}

#' Priors for transitions
#'
#' This function returns a keyed list of priors related to progression to
#' various health states. Called with no arguments, the default values are
#' returned. All parameters must be non-negative real numbers. 
#'
#' @param p_sym_if_inf A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of being
#'   symptomatic if infectious. 
#'
#'   Sources for default value: 
#'   
#'   \itemize{
#'     \item \insertRef{byambasuren_estimating_2020}{covidestim}
#'     \item \insertRef{mizumoto_estimating_2020}{covidestim}
#'     \item \insertRef{nishiura_estimation_2020}{covidestim}
#'   }
#'
#' @param p_sev_if_sym A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of
#'   severe if symptomatic. 
#'
#'   Sources for default value: 
#'
#'   \itemize{
#'     \item \insertRef{cdc_covid-19_response_team_severe_2020}{covidestim}
#'     \item \insertRef{verity_estimates_2020}{covidestim}
#'   }
#'
#' @param p_die_if_sev A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of dying
#'   if severely ill. 
#'
#'   Source for default value: \insertRef{cdc_covid-19_response_team_severe_2020}{covidestim}
#'
#' @param p_die_if_inf A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of dying
#'   if infected (e.g. the infection fatality rate). This prior represents a 
#'   national average value, which is adjusted for state and county-level 
#'   factors, and to reflect higher fatality rates early in the epidemic.
#'
#'   Source for default value: \insertRef{ODriscoll_Nature_2020}{covidestim}
#'   
#' @param ifr_decl_OR A two-element numeric vector containing \code{c(shape,
#'   rate)} parameters of a Gamma distribution modeling the elevated IFR in early
#'   2020, relative to the present. Default value represents an IFR 30% higher
#'   in March of 2020, with 95% interval 10-50% higher.
#'
#'   Source for default value: Expert opinion
#'
#' @return An S3 object of class \code{priors}
#' @examples
#' cfg <- covidestim(ndays = 50) + priors_transitions(p_sym_if_inf = c(0.5, 0.2))
#' @export
priors_transitions <- function(p_sym_if_inf = c(5.1430, 3.5360),    # a/b 
                               p_sev_if_sym = c(1.8854, 20.002),    # a/b
                               p_die_if_sev = c(28.239, 162.30),    # a/b
                               p_die_if_inf = c(15.915,3167.1),     # a/b
                               ifr_decl_OR  = c(9.1357, 29.339)     # a/b this is actually a gamma distribution!!
                               ) {                                 

  att(length(p_sym_if_inf) == 2)
  att(length(p_sev_if_sym) == 2)
  att(length(p_die_if_sev) == 2)
  att(length(p_die_if_inf) == 2)
  att(length(ifr_decl_OR) == 2)
  att(is_nonNegativeReal(p_sym_if_inf))
  att(is_nonNegativeReal(p_sev_if_sym))
  att(is_nonNegativeReal(p_die_if_sev))
  att(is_nonNegativeReal(p_die_if_inf))
  att(is_nonNegativeReal(ifr_decl_OR))
  
  build_priors(
    p_sym_if_inf,
    p_sev_if_sym,
    p_die_if_sev,
    p_die_if_inf,
    ifr_decl_OR,
    .postfix=c("_a", "_b"),
    .prefix="pri_"
  ) -> ps

  structure(ps, class='priors')
}

#' Fixed values on delay to progression
#'
#' This function returns a keyed list of values related to progression to
#' the infectious state. Called with no arguments, the default values are
#' returned. The following arguments can be passed to create different priors:
#'
#' @param inf_prg_delay A two-element numeric vector containing \code{c(shape, rate)}
#' parameters of a Gamma distribution modeling the time from infection to 
#' sypmtom onset. 
#'
#' Source for default value: \insertRef{lauer_incubation_2020}{covidestim}
#'
#' @param sym_prg_delay A two-element numeric vector containing \code{c(shape, rate)}
#' parameters of a Gamma distribution modeling the time from symptom onset to 
#' severe disease. 
#'
#' Source for default value: \insertRef{zhou_clinical_2020}{covidestim}
#'
#' @param sev_prg_delay A two-element numeric vector containing
#'   \code{c(shape, rate)} parameters of a Gamma distribution modeling the
#'   time from severe symptoms to death. 
#'   
#' Source for default value: \insertRef{linton_incubation_2020}{covidestim}
#' 
#' @param asy_rec_delay A two-element numeric vector containing 
#'    \code{c(shape,rate)} parameters of a Gamma distribution modeling the  
#'    time from infection to recovery without symptom development. 
#'    
#' @param pri_serial_i A two-element numeric vector containing 
#'    \code{c(shape,rate)} parameters of a Gamma distribution modeling the 
#'    serial interval, the average time between successive cases. 
#' 
#' @param infect_dist A two-element numeric vector containing 
#'    \code{c(shape,rate)} parameters of a Gamma distribution modeling the 
#'    changes of infectiousness, following initial infection.
#'    
#' @param seropos_dist A two-element numeric vector containing 
#'    \code{c(shape,rate)} parameters of a Gamma distribution modeling the 
#'    time to seroreversion, following initial infection. 
#' 
#' @return An S3 object of class \code{priors}
#' @examples
#' cfg <- covidestim(ndays = 50) + priors_progression(inf_prg_delay = c(4, 1))
#' @export
priors_progression <- function(inf_prg_delay = c(3.413, 0.6051), # shap/rate
                               sym_prg_delay = c(1.624, 0.2175), # shap/rate 
                               sev_prg_delay = c(2.061, 0.2277),  # shap/rate
                               asy_rec_delay = c(14   , 2     ),  # shap/rate 
                               pri_serial_i  = c(129.1, 22.25 ),  # shap/rate 
                               infect_dist   = c(8    , 1.241 ),  # shap/rate 
                               seropos_dist  = c(4.41 , 0.042 ) ) {   # shap/rate 

  att(length(inf_prg_delay) == 2)
  att(length(sym_prg_delay) == 2)
  att(length(sev_prg_delay) == 2)
  att(length(asy_rec_delay) == 2)
  att(length(pri_serial_i) == 2)
  att(length(infect_dist) == 2)
  att(length(seropos_dist) == 2)
  att(is_nonNegativeReal(inf_prg_delay))
  att(is_nonNegativeReal(sym_prg_delay))
  att(is_nonNegativeReal(sev_prg_delay))
  att(is_nonNegativeReal(asy_rec_delay))
  att(is_nonNegativeReal(pri_serial_i))
  att(is_nonNegativeReal(infect_dist))
  att(is_nonNegativeReal(seropos_dist))
  
  build_priors(
    inf_prg_delay,
    sym_prg_delay,
    sev_prg_delay,
    asy_rec_delay,
    pri_serial_i,
    infect_dist,
    seropos_dist,
    .postfix=c("_shap", "_rate")
  ) -> ps

  structure(ps, class="priors")
}

#' Priors on probability of diagnosis
#'
#' This function returns a keyed list of priors related to probability of
#' diagnosis. Called with no arguments, the default values are returned. 
#'
#' Boundary avoiding priors are used by default; a weakly informative prior
#' is used for the probability of diagnosis if severely ill. 
#'
#' @param rr_diag_asy_vs_sym A two-element numeric vector containing 
#'   \code{c(alpha, beta)} parameters of a Beta distribution modeling the rate 
#'   ratio of diagnosis among asymptomatic vs symptomatic, but not severe, 
#'   infections. 
#'   
#' @param rr_diag_sym_vs_asy A two-element numeric vector containing 
#'   \code{c(alpha, beta)} parameters of a Beta distribution modeling the rate 
#'   ratio of diagnosis at the symptomatic vs severe stage of infection. 
#'   
#' @param p_diag_if_sev A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of
#'   being diagnosed if severely ill. 
#'
#' @return An S3 object of class 'priors'
#' @examples
#' cfg <- covidestim(ndays = 50) + priors_diagnosis(p_diag_if_sym = c(2, 2))
#' @export
priors_diagnosis <- function(rr_diag_asy_vs_sym = c(2  ,18  ), # a/b
                             rr_diag_sym_vs_sev = c(2  , 2  ), # a/b
                             p_diag_if_sev      = c(5  , 2  )) {# a/b

  att(length(rr_diag_asy_vs_sym) == 2)
  att(length(rr_diag_sym_vs_sev) == 2)
  att(length(p_diag_if_sev) == 2)
  att(is_nonNegativeReal(rr_diag_asy_vs_sym))
  att(is_nonNegativeReal(rr_diag_sym_vs_sev))
  att(is_nonNegativeReal(p_diag_if_sev))

  build_priors(
    rr_diag_asy_vs_sym, 
    rr_diag_sym_vs_sev,
    p_diag_if_sev,
    .postfix=c("_a", "_b"),
    .prefix = "pri_"
  ) -> ps

  structure(ps, class="priors")
}

#' Fixed distributions on reporting delays
#'
#' This function returns a keyed list of mean values for the gamma distribution 
#' of delays related reporting of cases and deaths.  Called with no arguments, 
#' the default values are returned. 
#'
#' @param cas_rep_delay A two-element numeric vector containing the mean shape
#' and rate for a gamma distribution describing delay in case reporting. Shape
#' and rate are parameterized with a lognormal prior where the log of the input 
#' value is the median and the variance is 0.5. Based on average delay reported
#' in Santa Clara County, CA on 8 April 2020. 
#'
#' @param die_rep_delay A two-element numeric vector containing the mean shape
#' and rate for a gamma distribution describing delay in death reporting. Shape
#' and rate are parameterized with a lognormal prior where the log of the input 
#' value is the median and the variance is 0.5. Based on average delay reported
#' in Santa Clara County, CA on 8 April 2020. 
#'
#' @return An S3 object of class 'priors'
#' @examples
#' cfg <- covidestim(ndays = 50) + priors_reporting_delays(cas_rep_delay = c(0.5, 0.1))
#' @export
priors_reporting_delays <- function(cas_rep_delay = c(2.2,1),
                                    die_rep_delay = c(2.2,1)) {
  
  att(length(cas_rep_delay) == 2)
  att(length(die_rep_delay) == 2)
  att(is_nonNegativeReal(cas_rep_delay))
  att(is_nonNegativeReal(die_rep_delay))

  build_priors(
    cas_rep_delay,
    die_rep_delay,
    .postfix = c("_shap", "_rate")
  ) -> ps

  structure(ps, class = 'priors')
}

#' Priors on diagnostic delays
#'
#' This function returns a keyed list of priors related to delays in diagnosing
#' of cases and deaths.  Called with no arguments, the default values are 
#' returned. 
#'
#' Boundary avoiding priors are use by default. 
#' 
#' @param dx_delay_sym A two element vector containing the \code{c(alpha,beta)} 
#' parameters of a Beta distribution modeling a scaleing factor. Delay to 
#' diagnosis for symptomatic cases is modeled as the fraction of time in 
#' the symptomatic disease state, scaled by this factor. 
#'
#' @param dx_delay_sev A two element vector containing the \code{c(alpha,beta)} 
#' parameters of a Beta distribution modeling a scaleing factor. Delay to 
#' diagnosis for severe cases is modeled as the fraction of time in 
#' the severe disease state, scaled by this factor. 
#'
#' @return An S3 object of class 'priors'
#' @examples
#' cfg <- covidestim(ndays = 50) + priors_diagnosis_delays_scale(dx_delay_sym = c(0.5, 0.1))
#' @export
priors_diagnosis_delays_scale <- function(dx_delay_sym = c(2,2),
                                          dx_delay_sev = c(2,2)) {
  
  att(length(dx_delay_sym) == 2)
  att(length(dx_delay_sev) == 2)
  att(is_nonNegativeReal(dx_delay_sym))
  att(is_nonNegativeReal(dx_delay_sev))
  
  build_priors(
    dx_delay_sym,
    dx_delay_sev,
    .postfix = c("_a", "_b"),
    .prefix = "scale_"
  ) -> ps
  
  structure(ps, class = 'priors')
}

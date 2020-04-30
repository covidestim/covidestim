# Equations for converting a mean and variance to the alpha/beta parameters
# of a gamma distribution. Sourced from Reza's class slides.

#' Convert mean/variance stats to alpha/beta parameters of gamma distribution
#'
#' Using a method of moments approach, converts a sample with a specified
#' mean and variance to a parameterization of a gamma distribution
#'
#' @param mv A two-element numeric vector, where \code{[1]} represents mean
#'   and \code{[2]} represents variance.
#'
#' @return A number representing the alpha or beta parameter
pri_alpha <- function(mv) (mv[1]^2)/mv[2]

#' @rdname pri_alpha
pri_beta  <- function(mv) (mv[2]^2)/mv[1]

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
#' BUG: the following description needs information about what the defaults
#' were sourced from
#'
#' This function returns a keyed list of priors related to progression to
#' various health states. Called with no arguments, the default values are
#' returned. All parameters must be non-negative real numbers. Additionally,
#' the mean of the distribution described by \code{p_sym_if_inf} should be
#' \code{<=} than the mean of the
#' \code{\link[=priors_diagnosis]{p_diag_if_sym}}. Otherwise, a warning will be
#' triggered.
#'
#' @param p_sym_if_inf A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of being
#'   symptomatic if infectious
#' @param p_hos_if_sym A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of
#'   hospitalization if infectious
#' @param p_die_if_hos A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of dying
#'   if hospitalized
#'
#' @return An S3 object of class \code{priors}
#' @examples
#' cfg <- covidcast() + priors_transitions(p_sym_if_inf = c(0.5, 0.2))
#' @export
priors_transitions <- function(p_sym_if_inf = c(59, 41),      # a/b
                               p_hos_if_sym = c(31, 69),      # a/b
                               p_die_if_hos = c(1.25, 48.75)) { # a/b

  att(length(p_sym_if_inf) == 2)
  att(length(p_hos_if_sym) == 2)
  att(length(p_die_if_hos) == 2)
  att(is_nonNegativeReal(p_sym_if_inf))
  att(is_nonNegativeReal(p_hos_if_sym))
  att(is_nonNegativeReal(p_die_if_hos))

  build_priors(
    p_sym_if_inf,
    p_hos_if_sym,
    p_die_if_hos,
    .postfix=c("_a", "_b"),
    .prefix="pri_"
  ) -> ps

  structure(ps, class='priors')
}

#' Fixed values on delay to progression
#'
#' BUG: the following description needs information about what the defaults
#' were sourced from
#'
#' This function returns a keyed list of values related to progression to
#' the infectious state. Called with no arguments, the default values are
#' returned. The following arguments can be passed to create different priors:
#'
#' @param inf_prg_delay A two-element numeric vector containing \code{c(shape, scale)}
#' parameters of a Gamma distribution modeling [...]
#'
#' @param sym_prg_delay A two-element numeric vector containing \code{c(shape, scale)}
#' parameters of a Gamma distribution modeling [...] 
#'
#' @param hos_prg_delay A two-element numeric vector containing \code{c(shape, scale)}
#' parameters of a Gamma distribution modeling [...] 
#' 
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class \code{priors}
#' @examples
#' cfg <- covidcast() + priors_progression(inf_prg_delay = c(4, 1))
#' @export
priors_progression <- function(inf_prg_delay = c(5.202, 0.946), # shap/rate
                               sym_prg_delay = c(5.147, 0.468), # shap/rate 
                               hos_prg_delay = c(2.383, 0.27)) {# shap/rate

  att(length(inf_prg_delay) == 2)
  att(length(sym_prg_delay) == 2)
  att(length(hos_prg_delay) == 2)
  att(is_nonNegativeReal(inf_prg_delay))
  att(is_nonNegativeReal(sym_prg_delay))
  att(is_nonNegativeReal(hos_prg_delay))

  build_priors(
    inf_prg_delay,
    sym_prg_delay,
    hos_prg_delay,
    .postfix=c("_shap", "_rate")
  ) -> ps

  structure(ps, class="priors")
}

#' Fixed values on delay to recovered
#'
#' inf to res not well supported in data, guess:
#' assume mean 10 days, CI 5, 13
#' based on 5 days inf -> sym, 5 days sym -> res (CI 2,8)
#' adjust later?
#' BUG: ^ This needs to be integrated into the documentation!
#'
#' This function returns a keyed list of values related to recovery from
#' the infectious state. Called with no arguments, the default values are
#' returned. The following arguments can be passed to create different priors:
#'
#' @param inf_res_delay A two-element numeric vector containing \code{c(shape, scale)}
#' parameters of a Gamma distribution modeling [...]
#'
#' @param sym_res_delay A two-element numeric vector containing \code{c(shape, scale)}
#' parameters of a Gamma distribution modeling [...] 
#'
#' @param hos_res_delay A two-element numeric vector containing \code{c(shape, scale)}
#' parameters of a Gamma distribution modeling [...] 
#' 
#' @return An S3 object of class 'modelconfig'
#' @examples
#' cfg <- covidcast() + priors_recovery(inf_res_delay = c(20.12, 2.1))
#' @export
priors_recovery <- function(inf_res_delay = c(23.83, 2.383), # shap/rate
                            sym_res_delay = c(10.50, 2.099), # shap/rate
                            hos_res_delay = c(60.86, 3.567)) {# shap/rate

  att(length(inf_res_delay) == 2)
  att(length(sym_res_delay) == 2)
  att(length(hos_res_delay) == 2)
  att(is_nonNegativeReal(inf_res_delay))
  att(is_nonNegativeReal(sym_res_delay))
  att(is_nonNegativeReal(hos_res_delay))

  build_priors(
    inf_res_delay,
    sym_res_delay,
    hos_res_delay,
    .postfix=c("_shap", "_rate")
  ) -> ps

  structure(ps, class="priors")
}

#' Priors on reporting delay
#'
#' This function returns a keyed list of priors related to reporting delay.
#' Called with no arguments, the default values are returned. The following
#' arguments can be passed to create different priors:
#'
#' @param case_rep_delay A two-element numeric vector containing \code{c(shape,
#'   scale)} parameters of a Gamma distribution modeling [...]
#'
#' @param hos_rel_delay A two-element numeric vector containing \code{c(shape,
#'   scale)} parameters of a Gamma distribution modeling [...] 
#'
#' @param die_rep_delay A two-element numeric vector containing \code{c(shape,
#'   scale)} parameters of a Gamma distribution modeling [...] 
#'
#' All parameters should a \code{shape/rate < N_days}, where \code{N_days} is
#' the number of days of data being modeled, as passed to
#' \code{\link{covidcast()}}.
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' cfg <- covidcast() + priors_reporting_delay(pri_report_delay_shap = 6)
#' @export
priors_reporting_delay <- function(cas_rep_delay = c(1.73, 0.78), # shap/rate 
                                   hos_rep_delay = c(1.73, 0.78), # shap/rate 
                                   die_rep_delay = c(1.73, 0.78)) { # shap/rate

  att(length(cas_rep_delay) == 2)
  att(length(hos_rep_delay) == 2)
  att(length(die_rep_delay) == 2)
  att(is_nonNegativeReal(cas_rep_delay))
  att(is_nonNegativeReal(hos_rep_delay))
  att(is_nonNegativeReal(die_rep_delay))

  build_priors(
    cas_rep_delay,
    hos_rep_delay,
    die_rep_delay,
    .postfix=c("_shap", "_rate"),
    .prefix = "pri_"
  ) -> ps

  structure(ps, class="priors")
}

#' Priors on probability of diagnosis
#'
#' This function returns a keyed list of priors related to probability of
#' diagnosis.  Called with no arguments, the default values are returned. 
#'
#' Default values should be explained here, as well as constraints on custom
#' values. (Default values are assumed?)
#'
#' @param p_diag_if_inf A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of being
#'   diagnosed if infectious
#'
#' @param p_diag_if_sym A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of
#'   diagnosed if symptomatic. The mean of this distribution should be less
#'   the mean of the distribution parameterized by \code{p_diag_if_hos}.
#'   Otherwise, a warning will be displayed.
#'
#' @param p_diag_if_hos A two-element numeric vector containing \code{c(alpha,
#'   beta)} parameters of a Beta distribution modeling the probability of
#'   being diagnosed if hospitalized
#'
#' @return An S3 object of class 'priors'
#' @examples
#' cfg <- covidcast() + priors_diagnosis(p_diag_if_inf = c(0.5, 0.1))
#' @export
priors_diagnosis <- function(p_diag_if_inf = c(1, 99), # a/b
                             p_diag_if_sym = c(3, 7), # a/b
                             p_diag_if_hos = c(8.5, 1.5)) {# a/b

  att(length(p_diag_if_inf) == 2)
  att(length(p_diag_if_sym) == 2)
  att(length(p_diag_if_hos) == 2)
  att(is_nonNegativeReal(p_diag_if_inf))
  att(is_nonNegativeReal(p_diag_if_sym))
  att(is_nonNegativeReal(p_diag_if_hos))

  build_priors(
    p_diag_if_inf,
    p_diag_if_sym,
    p_diag_if_hos,
    .postfix=c("_a", "_b"),
    .prefix = "pri_"
  ) -> ps

  structure(ps, class="priors")
}

#' Priors on reporting delays
#'
#' This function returns a keyed list of priors related to delays in reporting
#' of cases, hospitalizations, and deaths.  Called with no arguments, the
#' default values are returned. 
#'
#' Default values should be explained here, as well as constraints on custom
#' values. (Default values are assumed?)
#'
#' @param cas_rep_delay Documentation about this prior
#'
#' @param hos_rep_delay Documentation about this prior
#'
#' @param die_rep_delay Documentation about this prior
#'
#' @return An S3 object of class 'priors'
#' @examples
#' cfg <- covidcast() + priors_reporting_delays_new(p_diag_if_inf = c(0.5, 0.1))
#' @export
priors_reporting_delays_new <- function(cas_rep_delay = c(3,1),
                                        hos_rep_delay = c(3,1),
                                        die_rep_delay = c(3,1)) {
  
  att(length(cas_rep_delay) == 2)
  att(length(hos_rep_delay) == 2)
  att(length(die_rep_delay) == 2)
  att(is_nonNegativeReal(cas_rep_delay))
  att(is_nonNegativeReal(hos_rep_delay))
  att(is_nonNegativeReal(die_rep_delay))

  build_priors(
    cas_rep_delay,
    hos_rep_delay,
    die_rep_delay,
    .postfix = c("_shp_a", "_shp_b")
  ) -> ps

  structure(ps, class = 'priors')
}

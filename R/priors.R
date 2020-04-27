# Equations for converting a mean and variance to the alpha/beta parameters
# of a gamma distribution. Sourced from Reza's class slides.
alpha_ <- function(mv) (mv[1]^2)/mv[2]
beta_  <- function(mv) (mv[2]^2)/mv[1]

# A list with arbitrary values, of class 'priors'
priors <- function(...) structure(list(...), class='priors')

print.priors <- function(ps, .tab = FALSE) {
'Priors:

' -> msg

  cat(msg)

  for (idx in names(ps)) {
    if (str_detect(idx, '^pri')) {
      idx_better <- str_replace(idx, '^pri_', '')
      idx_better <- str_replace(idx_better, '_a', '\t[alpha]')
      idx_better <- str_replace(idx_better, '_b', '\t[beta]')
      idx_better <- str_replace(idx_better, '_shap', '\t[shape]')
      idx_better <- str_replace(idx_better, '_rate', '\t[rate]')

      if (.tab)
        cat('\t')

      cat(glue("{idx_better}\t{ps[[idx]]}\n\n"))
    }
  }
  cat('\n')
}

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

#' @importFrom magrittr %>%
build_priors <- function(..., .postfix = c("_a", "_b")) {

  arg_names <- purrr::map_chr(rlang::enquos(...), rlang::quo_name)

  args_rekeyed <- stats::setNames(list(...), arg_names)

  purrr::lmap(
    args_rekeyed,
    ~rlang::dots_list(
      !! glue("pri_{names(.)}{.postfix[1]}") := .[[1]][1],
      !! glue("pri_{names(.)}{.postfix[2]}") := .[[1]][2],
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
#' returned. 
#'
#' @param p_sym_if_inf A two-element numeric vector containing mean, and
#'   variance of the probability of being symtomatic if infectious
#' @param p_hos_if_sym A two-element numeric vector containing mean, and
#'   variance of the probability of hospitalization if symptomatic
#' @param p_die_if_hos A two-element numeric vector containing mean, and
#'   variance of the probability of dying if hospitalized
#'
#' @return An S3 object of class \code{priors}
#' @examples
#' setup <- covidcast() + priors_transitions(p_sym_if_inf = c(0.5, 0.2))
#' @export
priors_transitions <- function(p_sym_if_inf = c(50, 50),   # a/b
                               p_hos_if_sym = c(30, 70),   # a/b
                               p_die_if_hos = c(2.5, 97.5)) { # a/b

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
    .postfix=c("_a", "_b")
  ) -> ps

  structure(ps, class='priors')
}

#' Priors on delay to progression
#'
#' BUG: the following description needs information about what the defaults
#' were sourced from
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
#' setup <- covidcast() + priors_progression(inf_prg_delay_shap = 3.5)
priors_progression <- function(inf_prg_delay = c(5.202, 0.946), # shap/rate
                               sym_prg_delay = c(5.147, 0.468), # shap/rate 
                               hos_prg_delay = c(9.164, 1.041)) {# shap/rate

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

#' Priors on delay to recovered
#'
#' inf to res not well supported in data, guess:
#' assume mean 10 days, CI 5, 13
#' based on 5 days inf -> sym, 5 days sym -> res (CI 2,8)
#' adjust later?
#' BUG: ^ This needs to be integrated into the documentation!
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
#' setup <- covidcast() + priors_recovery(inf_res_delay_shap = 20.12)
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
#' \code{pri_cas_rep_delay_shap}
#' \code{pri_cas_rep_delay_rate}
#' \code{pri_hos_rep_delay_shap}
#' \code{pri_hos_rep_delay_rate}
#' \code{pri_die_rep_delay_shap}
#' \code{pri_die_rep_delay_rate}
#'
#' Default values should be explained here, as well as constraints on custom
#' values. (Default values are assumed?)
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' setup <- covidcast() + priors_reporting_delay(pri_report_delay_shap = 6)
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
    .postfix=c("_shap", "_rate")
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
#' @param p_diag_if_inf A two-element numeric vector containing mean, and
#'   variance of the probability of being diagnosed if infectious
#' @param p_diag_if_sym A two-element numeric vector containing mean, and
#'   variance of the probability of being diagnosed if symptomatic
#' @param p_diag_if_hos A two-element numeric vector containing mean, and
#'   variance of the probability of being diagnosed if hospitalized
#'
#' @return An S3 object of class 'priors'
#' @examples
#' setup <- covidcast() + priors_diagnosis(p_diag_if_inf = c(0.5, 0.1))
#' @export
priors_diagnosis <- function(p_diag_if_inf = c(0.1, 9.9), # a/b
                             p_diag_if_sym = c(8.0, 2.0), # a/b
                             p_diag_if_hos = c(9.5, 0.5)) {# a/b

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
    .postfix=c("_a", "_b")
  ) -> ps

  structure(ps, class="priors")
}

priors_fixed <- function() {
  list(
    inf_prg_delay_shap_a = 4, 
    inf_prg_delay_shap_b = 1, 
    sym_prg_delay_shap_a = 4, 
    sym_prg_delay_shap_b = 1, 
    hos_prg_delay_shap_a = 4,
    hos_prg_delay_shap_b = 1,
    inf_res_delay_shap_a = 4,
    inf_res_delay_shap_b = 1,
    sym_res_delay_shap_a = 4,
    sym_res_delay_shap_b = 1,
    hos_res_delay_shap_a = 4,
    hos_res_delay_shap_b = 1,
    cas_rep_delay_shp_a  = 3,
    cas_rep_delay_shp_b  = 1.5, 
    hos_rep_delay_shp_a  = 3, 
    hos_rep_delay_shp_b  = 1.5,
    die_rep_delay_shp_a  = 3, 
    die_rep_delay_shp_b  = 1.5
  ) -> defaults

  structure(defaults, class = 'priors')
}

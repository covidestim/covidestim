#' Priors for transitions
#'
#' BUG: the following description needs information about what the defaults
#' were sourced from
#'
#' This function returns a keyed list of priors related to progression to
#' various health states. Called with no arguments, the default values are
#' returned. The following arguments can be passed to create different priors:
#'
#' \code{inf_prg_delay_shap}
#' \code{inf_prg_delay_rate}
#' \code{sym_prg_delay_shap}
#' \code{sym_prg_delay_rate}
#' \code{hos_prg_delay_shap}
#' \code{hos_prg_delay_rate}
#' 
#' Default values of these priors go here.
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' myData <- defaultConfig + priors_progression(inf_prg_delay_shap = 3.5)
#' @export
priors_transitions <- function(...) {
  args <- list(...)

  list(
    pri_p_sym_if_inf_a = 69,
    pri_p_sym_if_inf_b = 31,
    pri_p_hos_if_sym_a = 31,
    pri_p_hos_if_sym_b = 69,
    pri_p_die_if_hos_a = 3,
    pri_p_die_if_hos_b = 97
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
}

# A list with arbitrary values, of class 'priors'
priors <- function(...) structure(list(...), class='priors')

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
#' myData <- defaultConfig + priors_progression(inf_prg_delay_shap = 3.5)
#' @export
priors_progression <- function(...) {
  args <- list(...)

  list(
    inf_prg_delay_shap = 5.202,
    inf_prg_delay_rate = 0.946,
    sym_prg_delay_shap = 51.47,
    sym_prg_delay_rate = 4.679,
    hos_prg_delay_shap = 91.64,
    hos_prg_delay_rate = 10.41
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
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
#' myData <- defaultConfig + priors_recovery(inf_res_delay_shap = 20.12)
#' @export
priors_recovery <- function(...) {
  args <- list(...)

  list(
    inf_res_delay_shap = 23.83,
    inf_res_delay_rate = 2.383,
    sym_res_delay_shap = 10.50,
    sym_res_delay_rate = 2.099,
    hos_res_delay_shap = 60.86,
    hos_res_delay_rate = 3.567
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
}

#' Priors on reporting delay
#'
#' This function returns a keyed list of priors related to reporting delay.
#' Called with no arguments, the default values are returned. The following
#' arguments can be passed to create different priors:
#'
#' \code{pri_report_delay_shap}
#' \code{pri_report_delay_rate}
#'
#' Default values should be explained here, as well as constraints on custom
#' values. (Default values are assumed?)
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' myData <- defaultConfig + priors_reporting_delay(pri_report_delay_shap = 6)
#' @export
priors_reporting_delay <- function(...) {
  args <- list(...)

  list(
    pri_report_delay_shap = 7,
    pri_report_delay_rate = 0.9
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
}

#' Priors on probability of diagnosis
#'
#' This function returns a keyed list of priors related to probability of
#' diagnosis.  Called with no arguments, the default values are returned. The
#' following arguments can be passed to create different priors:
#'
#' \code{pri_p_diag_if_inf_a}
#' \code{pri_p_diag_if_inf_b}
#' \code{pri_p_diag_if_sym_a}
#' \code{pri_p_diag_if_sym_b}
#' \code{pri_p_diag_if_hos_a}
#' \code{pri_p_diag_if_hos_b}
#'
#' Default values should be explained here, as well as constraints on custom
#' values. (Default values are assumed?)
#'
#' @param ... A set of key-value pairs from those listed above. If no
#'   arguments are passed, default values will be used
#'
#' @return An S3 object of class 'modelconfig'
#' @examples
#' myData <- defaultConfig + priors_diagnosis(pri_p_diag_if_inf_a = 1.2)
#' @export
priors_diagnosis <- function(...) {
  args <- list(...)

  list(
    pri_p_diag_if_inf_a = 1,
    pri_p_diag_if_inf_b = 99,
    pri_p_diag_if_sym_a = 60,
    pri_p_diag_if_sym_b = 40,
    pri_p_diag_if_hos_a = 95,
    pri_p_diag_if_hos_b = 5
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
}


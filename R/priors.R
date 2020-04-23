# Equations for converting a mean and variance to the alpha/beta parameters
# of a gamma distribution. Sourced from Reza's class slides.
alpha_ <- function(mv) (mv[1]^2)/mv[2]
beta_  <- function(mv) (mv[2]^2)/mv[1]

#' Substitute a custom prior
#'
#' @param defaults The list of default priors
#' @param name A string. The prefix of the prior being modified. That is, 
#'   the name of the key of the prior, minus the '_a' and '_b' part.
#' @param mv_pair 2-element numeric. The pair of mean and variance that the
#'   user wants to substitute.
#'
#' @return The defaults, with the substitutes transfomed into alpha/beta
#'   parameters. If \code{mv_pair} was null, returns \code{defaults}
#'   unmodified.
subst_custom_prior <- function(defaults, name, mv_pair) {
  if (is.null(mv_pair))
    return(defaults)

  att(length(mv_pair) == 2)
  att(is.numeric(mv_pair))

  defaults[[paste0(name, "_a")]] <- alpha_(mv_pair)
  defaults[[paste0(name, "_b")]] <-  beta_(mv_pair)

  message(glue("Set custom prior for {name}:"))
  message(glue("mean={mv_pair[1]}, variance={mv_pair[2]}, alpha={alpha_(mv_pair)}, beta={beta_(mv_pair)}"))

  defaults
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
#' setup <- covidcast() + priors_progression(p_sym_if_inf = c(0.5, 0.2))
#' @export
priors_transitions <- function(p_sym_if_inf = NULL, p_hos_if_sym = NULL,
                               p_die_if_hos = NULL) {

  list(
    pri_p_sym_if_inf_a = 50,
    pri_p_sym_if_inf_b = 50,
    pri_p_hos_if_sym_a = 30,
    pri_p_hos_if_sym_b = 70,
    pri_p_die_if_hos_a = 2.5,
    pri_p_die_if_hos_b = 97.5
  ) -> defaults

  purrr::reduce2(
    c('pri_p_sym_if_inf',
      'pri_p_hos_if_sym',
      'pri_p_die_if_hos'),
      list(p_sym_if_inf,
           p_hos_if_sym,
           p_die_if_hos),
      subst_custom_prior,
      .init=list()) -> substitutions

  if (length(substitutions) > 0)
    att(all(names(substitutions) %in% names(defaults)))

  splice_class(defaults, substitutions, 'priors')
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
#' setup <- covidcast() + priors_progression(inf_prg_delay_shap = 3.5)
priors_progression <- function() {
  args <- list()

  list(
    pri_inf_prg_delay_shap = 5.202, 
    pri_nf_prg_delay_rate  = 0.946,
    pri_sym_prg_delay_shap = 5.147,
    pri_sym_prg_delay_rate = 0.468,
    pri_hos_prg_delay_shap = 9.164,
    pri_hos_prg_delay_rate = 1.041
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
#' setup <- covidcast() + priors_recovery(inf_res_delay_shap = 20.12)
priors_recovery <- function() {
  args <- list()

  list(
    pri_inf_res_delay_shap = 23.83,
    pri_inf_res_delay_rate = 2.383,
    pri_sym_res_delay_shap = 10.50,
    pri_sym_res_delay_rate = 2.099,
    pri_hos_res_delay_shap = 60.86,
    pri_hos_res_delay_rate = 3.567
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
priors_reporting_delay <- function() {
  args <- list()

  list(
    pri_cas_rep_delay_shap = 1.73,
    pri_cas_rep_delay_rate = 0.78,
    pri_hos_rep_delay_shap = 1.73, 
    pri_hos_rep_delay_rate = 0.78, 
    pri_die_rep_delay_shap = 1.73, 
    pri_die_rep_delay_rate = 0.78
  ) -> defaults

  att(all(names(args) %in% names(defaults)))

  splice_class(defaults, args, 'priors')
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
priors_diagnosis <- function(p_diag_if_inf = NULL, p_diag_if_sym = NULL,
                             p_diag_if_hos = NULL) {

  list(
    pri_p_diag_if_inf_a = 0.1,
    pri_p_diag_if_inf_b = 9.9,
    pri_p_diag_if_sym_a = 8.0,
    pri_p_diag_if_sym_b = 2.0,
    pri_p_diag_if_hos_a = 9.5,
    pri_p_diag_if_hos_b = 0.5
  ) -> defaults

  purrr::reduce2(
    c('pri_p_diag_if_inf',
      'pri_p_diag_if_sym',
      'pri_p_diag_if_hos'),
      list(p_diag_if_inf,
           p_diag_if_sym,
           p_diag_if_hos),
      subst_custom_prior,
      .init=list()) -> substitutions

  if (length(substitutions) > 0)
    att(all(names(substitutions) %in% names(defaults)))

  splice_class(defaults, substitutions, 'priors')
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

  splice_class(defaults, list(), 'priors')
}

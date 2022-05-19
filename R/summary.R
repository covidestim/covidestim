#' Summarize a Covidestim run
#'
#' Returns a \code{data.frame} summarizing a Covidestim model run. Note that if
#' \code{\link{runOptimizer.covidestim}} is used, all \code{*.(lo|hi)}
#' variables will be \code{NA}-valued, because BFGS does not generate
#' confidence intervals.
#'
#' @param ccr A \code{covidestim_result} object
#'
#' @param include.before A logical scalar. Include estimations that fall in the
#'   period before the first day of input data? (This period is of length
#'   \code{ndays_before} as passed to \code{covidestim}). If  \code{TRUE}, any
#'   elements of variables which do not have values for this "before" period
#'   will be represented as \code{NA}.
#'
#' @param index A logical scalar. If \code{TRUE}, will include a variable
#'   \code{index} in the output, with range \code{1:(ndays_before + ndays)}.
#'
#' @noMd
#' @return A \code{data.frame} with the following variables:
#'
#'   \itemize{
#'     \item \code{date}
#'
#'       Date as a \code{Date} vector.
#'
#'     \item \code{cases_fitted}, \code{cases_fitted_p2_5},
#'       \code{cases_fitted_p97_5}
#'
#'       The number of modeled case reports for date \code{date}. This is the
#'       model's expected number of filed case reports on date \code{date}.
#'       This number is computed downstream from the number of infections. It
#'       reflects underascertainment (not everyone infected is diagnosed),
#'       diagnosis delay (not everyone infected gets diagnosed) right
#'       away), and reporting delay (not all diagnoses are reported right
#'       away). (Cf. \strong{diagnoses} is the modeled number of cases
#'       (diagnosed infections) on date \code{date}, while
#'       \strong{cases_fitted} is the modeled number of \strong{reported} cases
#'       on date \code{date}).
#'
#'       \emph{Median and 95\% interval, ℝ}.
#
#'     \item \code{cumulative_immunoexposed},
#'       \code{cumulative_immunoexposed_p2_5},
#'       \code{cumulative_immunoexposed_p97_5}
#'
#'       Description of the parameter.
#'       
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{cumulative_infections}, \code{cumulative_infections_p2_5},
#'       \code{cumulative_infections_p97_5}
#'       
#'       The number of modeled cumulative infections (not cases, or diagnoses)
#'       that have occurred by the end of date \code{date}.
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'
#'     \item \code{deaths}, \code{deaths_p2_5}, \code{deaths_p97_5}
#'
#'       The number of modeled deaths for date \code{date}. The number of
#'       deaths estimated to have occurred on date \code{date} and does not
#'       account for reporting delays.
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'
#'     \item \code{deaths_diagnosed}, \code{deaths_diagnosed_p2_5},
#'       \code{deaths_diagnosed_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{deaths_fitted}, \code{deaths_fitted_p2_5},
#'       \code{deaths_fitted_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{deaths_fitted}, \code{deaths_fitted_p2_5},
#'       \code{deaths_fitted_p97_5}
#'
#'       The number of modeled death reports for a date \code{date}. This
#'       will always differ from the number of observed deaths for that same
#'       date, because \code{deaths_fitted} is approximating how many death
#'       reports should exist for that date.
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'
#'     \item
#'       \code{diagnoses}, \code{diagnoses_p2_5}, \code{diagnoses_p97_5}
#'
#'       The number of modeled diagnoses that occurred on date \code{date}.
#'       This is the sum of:
#'
#'       \itemize{
#'         \item New asymptomatic diagnoses on date \code{date}, \emph{plus:}
#'         \item New diagnoses of symptomatic, non-severe individuals on date
#'           \code{date}, \emph{plus:}
#'         \item New diagnoses of severe individuals on date \code{date}.
#'       }
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'
#'     \item \code{effective_protection_inf},
#'       \code{effective_protection_inf_p2_5},
#'       \code{effective_protection_inf_p97_5}
#'
#'       Description of the parameter.
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{effective_protection_inf_vax_boost},
#'       \code{effective_protection_inf_vax_boost_p2_5},
#'       \code{effective_protection_inf_vax_boost_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{effective_protection_inf_vax},
#'       \code{effective_protection_inf_vax_p2_5},
#'       \code{effective_protection_inf_vax_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{effective_protection_vax_boost},
#'       \code{effective_protection_vax_boost_p2_5},
#'       \code{effective_protection_vax_boost_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{effective_protection_vax},
#'       \code{effective_protection_vax_p2_5},
#'       \code{effective_protection_vax_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{first_inf}, \code{first_inf_p2_5}, \code{first_inf_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{fit_to_wastewater}, \code{fit_to_wastewater_p2_5},
#'       \code{fit_to_wastewater_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{infections}, \code{infections_p2_5}, \code{infections_p97_5}
#'
#'       The number of modeled infections that occurred on date \code{date}.
#'       This includes infections that may never cause symptoms, as well as
#'       infections which will never show up in case reports (will never be
#'       diagnosed). Being indexed by date-of-occurrence, reporting lag is
#'       absent from this outcome.  \emph{Median and 95\% interval, ℝ}.
#'
#'     \item \code{pop_susceptible}, \code{pop_susceptible_p2_5},
#'       \code{pop_susceptible_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{pop_susceptible_severe}, \code{pop_susceptible_severe_p2_5},
#'       \code{pop_susceptible_severe_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{Rt}, \code{Rt_p2_5}, \code{Rt_p97_5}
#'
#'       Estimate of the effective reproductive number (\eqn{R_t}).
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'     
#'     \item \code{sero_positive}, \code{sero_positive_p2_5}, 
#'       \code{sero_positive_p97_5}
#'
#'       The number of individuals in the population who are modeled as being
#'       seropositive. This is an unreliable outcome that we don't recommend
#'       using.
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'     
#'     \item \code{severe_fitted}, \code{severe_fitted_p2_5},
#'       \code{severe_fitted_p97_5}
#'
#'       Description of the parameter.
#'
#'       \emph{Median and 95\% interval, ℝ}
#'
#'     \item \code{severe}, \code{severe_p2_5}, \code{severe_p97_5}
#'
#'       The number of transitions into the "severe" health state on date
#'       \code{date}. The "severe" state is defined as disease that would merit
#'       hospitalization.  This outcome is not intended to model observational
#'       data detailing the number of COVID-positive hospital admissions or
#'       COVID-primary-cause hospital admissions.
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'
#'     \item \code{symptomatic_diagnosed}, \code{symptomatic_diagnosed_p2_5},
#'       \code{symptomatic_diagnosed_p97_5}
#'
#'       The number of modeled diagnoses of symptomatic individuals occurring
#'       on date \code{date}. The difference between this outcome and
#'       the \code{diagnoses} outcome is that \code{diagnoses} includes
#'       modeled diagnoses of asymptomatic individuals.
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'
#'     \item \code{symptomatic}, \code{symptomatic_p2_5}, \code{symptomatic_p97_5}
#'
#'       The number of modeled transitions of infected individuals into the
#'       infected, symptomatic health state on date \code{date}. This takes
#'       into account the probability of becoming symptomatic and the delay
#'       between infection and presentation of symptoms.
#'
#'       \emph{Median and 95\% interval, ℝ}.
#'
#'     \item \code{data_available}
#'
#'       Was there input data available for date \code{date}? This should be
#'       \code{TRUE}, except for the first month of estimates.
#'
#'       \emph{logical}.
#'
#' @export
#' @importFrom magrittr %>%
summary.covidestim_result <- function(ccr, include.before = TRUE, index = FALSE) {

  # Used for dealing with indices and dates
  nweeks        <- ccr$config$N_weeks
  nweeks_before <- ccr$config$N_weeks_before
  nweeks_total  <- nweeks_before + nweeks
  start_date    <- as.Date(ccr$config$first_date, origin = '1970-01-01')

  # Dates an index and converts it to a date
  toDate <- function(idx) start_date + lubridate::weeks(idx - 1 - nweeks_before)

  # Mappings between names in Stan and variable names in the output `df`
  c(
    "cases_fitted"                       = "cases_fitted",
    "cumulative_immunoexposed"           = "cumulative_immunoexposed",
    "cumulative_infections"              = "cumulative_infections",
    "deaths"                             = "deaths",
    "deaths_diagnosed"                   = "deaths_diagnosed",
    "deaths_fitted"                      = "deaths_fitted",
    "diagnoses"                          = "diagnoses",
    "effective_protection_inf"           = "effective_protection_inf",
    "effective_protection_inf_vax_boost" = "effective_protection_inf_vax_boost",
    "effective_protection_inf_vax"       = "effective_protection_vax",
    "effective_protection_vax_boost"     = "effective_protection_vax_boost",
    "effective_protection_vax"           = "effective_protection_vax",
    "first_inf"                          = "first_infections",
    "fit_to_wastewater"                  = "fit_to_wastewater",
    "infections"                         = "infections",
    "pop_susceptible"                    = "pop_susceptible",
    "pop_susceptible_severe"             = "pop_susceptible_severe",
    "Rt"                                 = "Rt",
    "sero_positive"                      = "sero_positive",
    "severe_fitted"                      = "severe_fitted",
    "severe"                             = "severe",
    "symptomatic_diagnosed"              = "symptomatic_diagnosed",
    "symptomatic"                        = "symptomatic"
  ) -> params

  
  if ("optimizer" %in% ccr$flags)
    return(summaryOptimizer(ccr, toDate, params, start_date))

  # Used for renaming quantiles output by Stan
  quantile_names <- c("2.5%" = "_p2_5", "50%" = "", "97.5%" = "_p97_5")

  # This creates the list of `pars` that gets passed to `rstan::summary`
  # by enumerating all combs of varnames and indices
  purrr::cross2(names(params), as.character(1:ndays_total)) %>% 
  purrr::map_chr(
    function(item)
      glue("{par}[{idx}]", par = item[[1]], idx = item[[2]])
  ) -> pars

  rstan::summary(
    ccr$result,
    pars = pars,
    probs = c(0.025, 0.5, 0.975)
  )$summary %>% # Get rid of the per-chain summaries by indexing into `$summary`
    as.data.frame %>%
    tibble::as_tibble(rownames = "parname") -> melted

  # These are the variables that are going to be selected from the melted
  # representation created above
  vars_of_interest <- c("variable", "date", names(quantile_names))

  # Join the melted representation to the array indices that have been split
  # into their name:idx components
  stan_extracts <- dplyr::left_join(
    melted,
    split_array_indexing(melted$parname),
    by = "parname"
  ) %>%
    # Reformat the dates, and rename some of the variable names
    dplyr::mutate(date = toDate(index), variable = params[variable]) %>%
    # Eliminate things like R-hat that we don't care about right now
    dplyr::select_at(vars_of_interest) %>%
    # Melt things more to get down to three columns
    tidyr::gather(names(quantile_names), key = "quantile", value = "value") %>%
    # Create the finalized names for the quantiles, and delete the now-unneeded
    # quantile variable
    dplyr::mutate(
      variable = paste0(variable, quantile_names[quantile]),
      quantile = NULL,
      
    ) %>%
    # Cast everything back out
    tidyr::spread(key = "variable", value = "value") %>%
    dplyr::mutate(data_available = date >= start_date)

  d <- stan_extracts


  if (index == TRUE)
    d <- dplyr::bind_cols(d, list(index = 1:ndays_total))

  # Get rid of 'before period' data if not requested
  if (include.before == FALSE)
    d <- dplyr::filter(d, date >= start_date)

  d
}

split_array_indexing <- function(elnames) {

  # Matches a valid indexed variable in Stan, capturing it into two groups
  regex <- '^([A-Za-z_][A-Za-z_0-9]+)\\[([0-9]+)\\]$'

  # Match everything into a data.frame
  captured <- stringr::str_match(elnames, regex)
  captured <- as.data.frame(captured, stringsAsFactors = FALSE)

  # Rename cols, coerce the index to a number
  colnames(captured) <- c("parname", "variable", "index")
  captured$index <- as.numeric(captured$index)

  captured
}

#' Summarize internal parameters of a Covidestim run
#'
#' Returns a \code{data.frame} with a information on the posterior of a few
#' key Stan parameters.
#'
#' @param ccr A \code{covidestim_result} object
#'
#' @return A \code{data.frame} with the following variables:
#'
#'   \itemize{
#'     \item \code{par}: Name of the Stan parameter
#'
#'     \item \code{2.5\%}, \code{50\%}, \code{97.5\%}: Quantiles for the posterior
#'       distribution of \code{par}.
#'   }
#'
#' @export
#' @importFrom magrittr %>%
summaryEpi <- function(ccr) {
  
  c(
    "p_sym_if_inf",
    "p_sev_if_sym",
    "p_die_if_sev",
    "p_die_if_sym",
    "p_diag_if_sym",
    "p_diag_if_sev"
  ) -> pars_of_interest

  # Used for renaming quantiles output by Stan
  quantile_names <- c("2.5%" = "_p2_5", "50%" = "median", "97.5%" = "_p97_5")

  rstan::summary(
    ccr$result,
    pars = pars_of_interest,
    probs = c(0.025, 0.5, 0.975)
  )$summary %>% # Get rid of the per-chain summaries by indexing into `$summary`
    as.data.frame %>%
    tibble::as_tibble(rownames = "par") -> melted

  # These are the variables that are going to be selected from the melted
  # representation created above
  vars_of_interest <- c("par", names(quantile_names))

  result <- dplyr::select_at(melted, vars_of_interest)

  result
}

summaryOptimizer <- function(ccr, toDate, params, start_date) {

  # Optimizer results don't have confidence intervals - however, to maintain
  # parity with the "default" version of the summary function (where the
  # sampler produces the results that are being processed), we need to have
  # NA-valued confidence intervals. This also makes sharing a DB table with
  # sampler-generated results much easier. The next few lines create the
  # neccessary "*_p2_5" and "*_p97_5" NA-valued columns
  params.lo <- paste0(params, "_p2_5")
  params.hi <- paste0(params, "_p97_5")

  nullParams <- as.list(rep(NA, 2*length(params))) %>%
    stats::setNames(c(params.lo, params.hi))

  # ccr$result is a list keyed on Stan param name, each value is a numeric
  # vector of equal length.
  tibble::as_tibble(ccr$result) %>%
    dplyr::mutate_all(as.double) %>%
    # Convert to result naming-scheme provided by `summary()`
    dplyr::rename_with(~params[.]) %>%
    # Add empty conf. interval variables
    dplyr::bind_cols(nullParams) %>%
    # Sort result columns in alphabetical order to maintain parity with
    # original `summary` implementation
    dplyr::select(order(colnames(.))) %>%
    # Add dates and data_available column
    dplyr::mutate(
      date = toDate(1:dplyr::n()),
      data_available = date >= start_date
    ) %>%
    # Maintain parity with order in original `summary` implementation
    dplyr::select(date, everything(), data_available)
}

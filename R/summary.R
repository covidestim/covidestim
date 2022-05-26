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
#'   \code{nweeks_before} as passed to \code{covidestim}). If  \code{TRUE}, any
#'   elements of variables which do not have values for this "before" period
#'   will be represented as \code{NA}.
#'
#' @param index A logical scalar. If \code{TRUE}, will include a variable
#'   \code{index} in the output, with range \code{1:(nweeks_before + nweeks)}.
#'
#' @noMd
#' @return A \code{data.frame} with the following variables:
#'
#'   \itemize{
#'     \item \bold{\code{date}}
#'
#'       Date as a \code{Date} vector.
#'
#'     \item \bold{\code{deaths}}, \code{deaths} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{deaths_of_diagnosed}}, \code{deaths_of_diagnosed} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{diagnoses}}, \code{diagnoses} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{diagnoses_of_symptomatic}}, \code{diagnoses_of_symptomatic} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{effective_protection_inf_pvl}}, \code{effective_protection_inf_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{effective_protection_inf_vax_boost_pvl}}, \code{effective_protection_inf_vax_boost_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{effective_protection_inf_vax_pvl}}, \code{effective_protection_inf_vax_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{effective_protection_vax_boost_pvl}}, \code{effective_protection_vax_boost_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{effective_protection_vax_pvl}}, \code{effective_protection_vax_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{fitted_cases}}, \code{fitted_cases} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{fitted_deaths}}, \code{fitted_deaths} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{fitted_hospitalizations}}, \code{fitted_hospitalizations} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{fitted_wastewater_pvl}}, \code{fitted_wastewater_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{immunoexposed_cumulative}}, \code{immunoexposed_cumulative} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{infections_cumulative}}, \code{infections_cumulative} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{infections}}, \code{infections} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{infections_premiere}}, \code{infections_premiere} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{r_t}}, \code{r_t} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{seropositive_pvl}}, \code{seropositive_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{severe}}, \code{severe} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{susceptible_pvl}}, \code{susceptible_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{susceptible_severe_pvl}}, \code{susceptible_severe_pvl} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{symptomatic}}, \code{symptomatic} + (\code{p2_5}, \code{p25}, \code{p75}, \code{p97_5})
#'
#'       Description of the parameter.
#'
#'       \emph{Median, 2.5\% interval, 25\% interval, 75\% interval, 97.5\%
#'         interval, RRR}.
#'
#'     \item \bold{\code{data_available}}
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
    "deaths"                                 = "deaths",
    "deaths_of_diagnosed"                    = "deaths_of_diagnosed",
    "diagnoses"                              = "diagnoses",
    "diagnoses_of_symptomatic"               = "diagnoses_of_symptomatic",
    "effective_protection_inf_pvl"           = "effective_protection_inf_pvl",
    "effective_protection_inf_vax_boost_pvl" = "effective_protection_inf_vax_boost_pvl",
    "effective_protection_inf_vax_pvl"       = "effective_protection_inf_vax_pvl",
    "effective_protection_vax_boost_pvl"     = "effective_protection_vax_boost_pvl",
    "effective_protection_vax_pvl"           = "effective_protection_vax_pvl",
    "fitted_cases"                           = "fitted_cases",
    "fitted_deaths"                          = "fitted_deaths",
    "fitted_hospitalizations"                = "fitted_hospitalizations",
    "fitted_wastewater_pvl"                  = "fitted_wastewater_pvl",
    "immunoexposed_cumulative"               = "immunoexposed_cumulative",
    "infections_cumulative"                  = "infections_cumulative",
    "infections"                             = "infections",
    "infections_premiere"                    = "infections_premiere",
    "r_t"                                    = "r_t",
    "seropositive_pvl"                       = "seropositive_pvl",
    "severe"                                 = "severe",
    "susceptible_pvl"                        = "susceptible_pvl",
    "susceptible_severe_pvl"                 = "susceptible_severe_pvl",
    "symptomatic"                            = "symptomatic"
  ) -> params

  
  if ("optimizer" %in% ccr$flags)
    return(summaryOptimizer(ccr, toDate, params, start_date))

  # Used for renaming quantiles output by Stan
  quantile_names <- c(
    "2.5%"  = "_p2_5",
    "25%"   = "_p25",
    "50%"   = "",
    "75%"   = "_p75",
    "97.5%" = "_p97_5"
  )

  # This creates the list of `pars` that gets passed to `rstan::summary`
  # by enumerating all combs of varnames and indices
  purrr::cross2(names(params), as.character(1:nweeks_total)) %>% 
  purrr::map_chr(
    function(item)
      glue("{par}[{idx}]", par = item[[1]], idx = item[[2]])
  ) -> pars

  rstan::summary(
    ccr$result,
    pars = pars,
    probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
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
      variable = paste0(variable, quantile_names[I(quantile)]),
      quantile = NULL
    ) %>%
    # Cast everything back out
    tidyr::spread(key = "variable", value = "value") %>%
    dplyr::mutate(data_available = date >= start_date)

  d <- stan_extracts


  if (index == TRUE)
    d <- dplyr::bind_cols(d, list(index = 1:nweeks_total))

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
  quantile_names <- c(
    "2.5%"  = "_p2_5",
    "25%"   = "_p25",
    "50%"   = "",
    "75%"   = "_p75",
    "97.5%" = "_p97_5"
  )

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
  params.2_5  <- paste0(params, "_p2_5")
  params.25   <- paste0(params, "_p25")
  params.75   <- paste0(params, "_p75")
  params.97_5 <- paste0(params, "_p97_5")

  nullParams <- as.list(rep(NA, 4*length(params))) %>%
    stats::setNames(c(params.2_5, params.25, params.75, params.97_5))

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

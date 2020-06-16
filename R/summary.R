#' Summarize a Covidestim run
#'
#' @param ccr A \code{covidestim_result} object
#'
#' @param include.before A logical scalar. Include estimations that fall in the
#'   period before the first day of input data? (This period is of length
#'   \code{N_days_before} as passed to \code{covidestim}). If  \code{TRUE}, any
#'   elements of variables which do not have values for this "before" period
#'   will be represented as \code{NA}.
#'
#' @param include.RtEstim A logical scalar. Include Rt estimates in output?
#'   Will add ~1 minute to runtime.
#'
#' @param index A logical scalar. If \code{TRUE}, will include a variable
#'   \code{index} in the output, with range \code{1:(N_days_before + N_days)}.
#'
#' @return A \code{data.frame} with the following variables:
#'
#'   \itemize{
#'     \item \code{date}
#'     \item \code{cases.fitted}, \code{cases.fitted.lo},
#'       \code{cases.fitted.hi}: Median and 95\% HDI around fitted
#'       case data.
#'
#'     \item \code{cum.incidence}, \code{cum.incidence.lo},
#'       \code{cum.incidence.hi}: Median and 95\% HDI around
#'       cumulative incidence data.
#'
#'     \item \code{deaths.fitted}, \code{deaths.fitted.lo},
#'       \code{deaths.fitted.hi}: Median and 95\% HDI around fitted
#'       deaths data.
#'
#'     \item \code{deaths}, \code{deaths.lo}, \code{deaths.hi}: Median and 95\%
#'     HDI around estimated deaths/day.
#'
#'     \item \code{infections}, \code{infections.lo}, \code{infections.hi}:
#'       Median and 95\% HDI around estimated infections/day.
#'
#'     \item \code{severe}, \code{severe.lo}, \code{severe.hi}: Median and 95\%
#'        HDI around "severe" estimate.
#'
#'     \item \code{symptomatic}, \code{symptomatic.lo}, \code{symptomatic.hi}:
#'     Median and 95\% HDI around estimate of quantity of individuals
#'     symtomatic on a particular day.
#'
#'     \item \code{data.available}: \code{TRUE/FALSE} for whether input data
#'       was available on that particular day.
#'
#'     \item \code{Rt}, \code{Rt.lo}, \code{Rt.hi}: Estimate of Rt, with 95\%
#'       HDI. See \code{\link{RtEst}}.
#'
#'     \item \code{NaiveRt}, \code{NaiveRt.lo}, \code{NaiveRt.hi}: Naive
#'       estimate of Rt, with 95% HDI. See \code{\link{RtNaiveEstim}}.
#'
#'     \item \code{index}
#'   }
#'
#' @export
#' @importFrom magrittr %>%
summary.covidestim_result <- function(ccr, include.before = TRUE, 
                                     include.RtEstim = TRUE, index = FALSE) {

  # Used for dealing with indices and dates
  ndays        <- ccr$config$N_days
  ndays_before <- ccr$config$N_days_before
  ndays_total  <- ndays_before + ndays
  start_date   <- as.Date(ccr$config$first_date, origin = '1970-01-01')

  # Mappings between names in Stan and variable names in the output `df`
  c(
    "new_inf"              = "infections",
    "occur_cas"            = "cases.fitted",
    "occur_die"            = "deaths.fitted",
    "cumulative_incidence" = "cum.incidence",
    "new_sym"              = "symptomatic",
    "new_sev"              = "severe",
    "new_die"              = "deaths"
  ) -> params

  # Used for renaming quantiles output by Stan
  quantile_names <- c("50%" = "")

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
    probs = 0.5
  )$summary %>% # Exclude per-chain summaries by indexing into `$summary`
    as.data.frame %>%
    tibble::as_tibble(rownames = "parname") -> melted

  # These are the variables that are going to be selected from the melted
  # representation created above
  vars_of_interest <- c("variable", "date", names(quantile_names))

  # Dates an index and converts it to a date
  toDate <- function(idx) start_date + lubridate::days(idx - 1 - ndays_before)

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
    dplyr::mutate(data.available = date >= start_date)

  d <- stan_extracts

  # Join in Rt estimates if requested
  if (include.RtEstim) {
    Rt   <- RtEst(ccr)
    Rt_n <- RtNaiveEstim(ccr)
    d <- dplyr::left_join(d, Rt, by = "date")
    d <- dplyr::left_join(d, Rt_n, by = "date")
  }

  if (index == TRUE)
    d <- dplyr::bind_cols(d, list(index = 1:ndays_total))

  # Calculate HDI's for all parameters and name them appropriately
  purrr::map2_dfc(
    names(params), params,
    function(stan_name, df_name_prefix) {

      # Result of this is a list with keys $lo and $hi
      interval <- hdi_for_ccr_indexed_par(ccr, stan_name)

      # Compose the final names for these HDIs
      df_names <- paste0(df_name_prefix, c(".lo", ".hi"))
      
      rlang::dots_list(!!df_names[1] := interval$lo,
                       !!df_names[2] := interval$hi)
    }
  ) %>% as.list -> hdis

  d <- dplyr::bind_cols(d, hdis)

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

# 'x' is a numeric vector
hdi <- function(x) {
  # order samples
  x <- x[order(x)]

  # calc samples in XX% interval. lets say XX = 95
  N <- length(x)
  n <- N*0.95

  # find index for min
  xx <- NULL
  for(i in 1:(N-n)) xx[i] <- diff(x[c(0,n)+i])

  xx_min <- which(xx==min(xx))

  # calc intervals
  hd_interval <- x[c(0,n)+xx_min]             

  list(lo = hd_interval[1], hi = hd_interval[2])
}

# covidestim_result, string -> list(lo, hi) of HDIs for a particular Stan param
hdi_for_ccr_indexed_par <- function(ccr, parname) {

  att(parname %in% names(ccr$extracted))

  samples <- ccr$extracted[[parname]]

  ndays        <- ccr$config$N_days
  ndays_before <- ccr$config$N_days_before
  ndays_total  <- ndays_before + ndays

  start_date   <- as.Date(ccr$config$first_date, origin = '1970-01-01')

  purrr::map_dfr(1:ndays_total, ~hdi(samples[, .])) %>% as.list
}

# covidestim_result, string -> list(lo, hi) of HDIs for Stan params that are
# NOT indexed
hdi_for_ccr_par <- function(ccr, parname) {

  att(parname %in% names(ccr$extracted))

  samples <- ccr$extracted[[parname]]

  hdi(samples)
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
  quantile_names <- c("2.5%" = ".lo", "50%" = "median", "97.5%" = ".hi")

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


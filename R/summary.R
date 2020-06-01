#' @export
summary <- function(ccr, include.before = TRUE, 
                    include.RtEstim = TRUE, iter = 500,
                    index = FALSE) UseMethod('summary')

#' @export
#' @importFrom magrittr %>%
summary.covidcast_result <- function(ccr, include.before = TRUE, 
                                     include.RtEstim = TRUE, iter = 500,
                                     index = FALSE) {

  # Used for dealing with indices and dates
  ndays        <- ccr$config$N_days
  ndays_before <- ccr$config$N_days_before
  ndays_total  <- ndays_before + ndays
  start_date   <- as.Date(ccr$config$first_date, origin = '1970-01-01')

  # Mappings between names in Stan and variable names in the output `df`
  c(
    "new_inf" = "infections.est",
    "occur_cas" = "cases.fitted",
    "occur_die" = "deaths.fitted",
    "cumulative_incidence" = "cum.incidence.est"
  ) -> params

  # Used for renaming quantiles output by Stan
  quantile_names <- c("2.5%" = ".lo", "50%" = "", "97.5%" = ".hi")

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
    Rt   <- RtEst(ccr, graph = FALSE)
    Rt_n <- RtNaiveEstim(ccr, graph = FALSE)
    d <- dplyr::left_join(d, Rt, by = "date")
    d <- dplyr::left_join(d, Rt_n, by = "date")
  }

  # Get rid of 'before period' data if not requested
  if (include.before == FALSE)
    d <- dplyr::filter(d, date >= start_date)

  # Column-bind high-density interval CI's for cumulative incidence
  # First, remove the traditionally-calculated lo and hi estimates
  d <- dplyr::select(d, -cum.incidence.est.lo, -cum.incidence.est.hi)
  d <- dplyr::bind_cols(d, cum_inc_hdi(ccr))

  if (index == TRUE)
    d <- dplyr::bind_cols(d, list(index = 1:ndays_total))

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

cum_inc_hdi <- function(ccr) {
  samples <- ccr$extracted$cumulative_incidence

  ndays        <- ccr$config$N_days
  ndays_before <- ccr$config$N_days_before
  ndays_total  <- ndays_before + ndays

  start_date   <- as.Date(ccr$config$first_date, origin = '1970-01-01')

  purrr::map_dfr(1:ndays_total, ~hdi(samples[, .])) %>%
    dplyr::rename(cum.incidence.est.lo = lo,
                  cum.incidence.est.hi = hi)
}

#' @export
RtEst <- function(...) UseMethod('RtEst')

#' Calculate Effective Reproduction Number from Covidestim output
#'
#' Creates a \code{data.frame} with Rt estimates over time and 95\% CIs.
#'
#' This code makes use of the \code{EpiEstim} package, described here:
#'
#' \insertRef{cori_new_2013}{covidestim}
#'
#' The method used in \code{EpiEstim} is \code{"parametric_si"}, and a 5-day
#' sliding window is used by default. \strong{Warning messages specifying that 
#' you are \code{"estimating R too early"} are expected and should be ignored}.
#'
#' @param ccr The result of calling \code{\link{run}}. An object of class
#'   \code{covidestim_result}.
#'   
#' @param window.length The size of the moving window over which Rt should be 
#' estimated; defaults to 5 days. This number must be odd, so that the
#' interval's true center always falls on a complete day. It also must be less
#' than half the number of data days, i.e., \code{window.length < ndays/2}.
#' 
#' @param mean.si The mean serial interval to use in the Rt estimate; defaults
#' to 4.7 days. 
#' 
#' @param std.si The standard deviation of the serial interval; defaults to 
#' 2.9 days.  
#'
#' @return A \code{data.frame} with columns \code{date}, \code{Rt},
#'   \code{Rt.lo}, \code{Rt.hi}, where \code{date} is a vector of
#'   \code{\link[base]{Date}} objects. Or, for \code{RtNaiveEstim}, columns
#'   \code{date}, \code{NaiveRt}, \code{NaiveRt.lo}, \code{NaiveRt.hi}
#'
#' @export
RtEst.covidestim_result <- function(ccr, window.length = 5, mean.si = 4.7,
                                    std.si = 2.9) {

  # Calculating Rt for the entire period, not just the portion of the period
  # for which we have data
  n_days        <- ccr$config$N_days
  n_days_before <- ccr$config$N_days_before
  n_days_tot    <- n_days_before + n_days

  att(window.length %% 2 == 1)
  att(window.length < n_days / 2)

  # Grab the estimated number of infections/day for the whole period and
  # round it to an integer value
  new_inf    <- ccr$extracted$new_inf %>% round

  # 'first_data_date' is the first date of the 'N_days' period
  # 'first_date' is the first day of the 'N_days_before' period
  first_data_date <- ccr$config$first_date %>% as.Date(origin = '1970-01-01')
  first_date      <- first_data_date - lubridate::days(n_days_before)

  EpiEstim::make_config(
    list(
      # EpiEstim requires that the first day considered is not the first day
      # of data, hence the '2'
      t_start = 2:(n_days_tot - window.length + 1),
      t_end   = (2 + window.length - 1):n_days_tot,
      mean_si = mean.si, 
      std_si  = std.si
    )
  ) -> RtEstimConfig

  # This calculates an Rt+bounds for a single time series of incidence data
  runEpiEstim <- function(new_inf, iteration) {

    EpiEstim::estimate_R(
      new_inf,
      method = "parametric_si",
      config = RtEstimConfig
    )$R -> R_est

    # '+1': The first day of our first day is the second day
    # '+2': The middle day of that first 5-day window is two days away from the
    #       first day of that first 5-day window
    first_day_centered <- first_date + lubridate::days(1 + 2)
    last_day_centered  <- first_day_centered + lubridate::days(nrow(R_est) - 1)

    list(
      date   = seq(first_day_centered, last_day_centered, by = '1 day'),
      median = R_est$`Median(R)`
    )
  }

  # Run 'runEpiEstim' for each row of the 'new_inf' matrix, where a row
  # represents one iteration.
  EpiEstim.result <-
    purrr::imap_dfr(1:nrow(new_inf), ~runEpiEstim(new_inf[.x, ], .y))

  # Summarize all of these iterations by taking the median of 'Rt', and
  # calculating all quantiles based off of that
  dplyr::summarize(
    dplyr::group_by(EpiEstim.result, date),
    Rt    = median(I(median)), # The 'I()' call preserves meaning of 'median'
    Rt.lo = quantile(I(median), 0.025),
    Rt.hi = quantile(I(median), 0.975)
  )
}

#' @export
RtNaiveEstim <- function(...) UseMethod('RtNaiveEstim')

#' @rdname RtEst.covidestim_result
#' @export
RtNaiveEstim.covidcast_result <- function(cc, window.length = 5, mean.si = 4.7,
                                          std.si = 2.9) {

  day_start <- 2
  day_end   <- day_start + window.length - 1
  n_days    <- cc$config$N_days
  ndb       <- cc$config$N_days_before

  att(window.length %% 2 == 1)
  att(window.length < n_days / 2)

  first_date <- as.Date(cc$config$first_date, origin = '1970-01-01')

  # The case data used in this model run
  cases <- cc$config$obs_cas

  EpiEstim::make_config(
    list(
      t_start = seq(day_start, n_days - window.length + 1), 
      t_end   = seq(day_end,   n_days),
      mean_si = mean.si, 
      std_si  = std.si
    )
  ) -> config

  EpiEstim::estimate_R(
    cases,
    method = "parametric_si",
    config = config
  ) -> Rt_Est

  result <- tibble::as_tibble(Rt_Est[["R"]])
  
  # Pull quantiles directly from the result of EpiEstim::estimate_R()
  dplyr::transmute(
    result,
    date = seq(first_date + lubridate::days(1+2),
               first_date + lubridate::days(n_days - 1 - 2),
               by = '1 day'),
    NaiveRt = `Median(R)`,
    NaiveRt.lo = `Quantile.0.025(R)`,
    NaiveRt.hi = `Quantile.0.975(R)`
  )
}


#' @export
#' @rdname RtEst.covidcast_result
RtEst <- function(...) UseMethod('RtEst')

#' Calculate the Effective Reprouction Number from CovidCast output
#'
#' Creates a figure with Rt estimates over time and 95% CI and undelrying df. 
#'
#' @param ccr The result of calling \code{\link{run}}. An object of class
#'   \code{covidcast_result}.
#'   
#' @param sample_fraction The fraction of iterations to sample for the Rt
#' calualation; defaults to 2/3 samples after warm-up.  
#' 
#' @param window.length The size of the moving window over which Rt should be 
#' estimated; defaults to 5 days. This number must be odd. 
#' 
#' @param mean.si The mean serial interval to use in the Rt estimate; defaults
#' to 4.7 days. 
#' 
#' @param std.si The standard deviation of the serial interval; defaults to 
#' 2.9 days.  
#'
#' @param pdf.name This will render a pdf to the workingvdirectory with name
#'  \code{covidcast_Rt.pdf}.
#'
#' @return an df with estimates
#'
#' @export
validate_Rt_input <- function(d,e) {
  pvec <- purrr::partial(paste, ...=, collapse = ', ')
  att(
    d %% 2 == 1,
    msg="'window.length' must be an odd number"
  )
  att(
    0 < e && e <= 1,
    msg="'sample_fraction' must be greater than 0 and less than or eqaul to 1"
  )
}

#' @export
RtEst.covidcast_result <- function(ccr, window.length = 5, mean.si = 4.7,
                                   std.si = 2.9, graph = TRUE) {

  # Calculating Rt for the entire period, not just the portion of the period
  # for which we have data
  n_days        <- ccr$config$N_days
  n_days_before <- ccr$config$N_days_before
  n_days_tot    <- n_days_before + n_days

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

    # message(glue::glue("Iteration {iteration} beginning"))

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
      date = seq(first_day_centered, last_day_centered, by = '1 day'),
      mean = R_est$`Mean(R)`,
      lo   = R_est$`Quantile.0.025(R)`,
      hi   = R_est$`Quantile.0.975(R)`
    )
  }

  # Run 'runEpiEstim' for each row of the 'new_inf' matrix, where a row
  # represents one iteration.
  EpiEstim.result <-
    purrr::imap_dfr(1:nrow(new_inf), ~runEpiEstim(new_inf[.x, ], .y))

  # Summarize all of these iterations by taking the mean of 'Rt', the 2.5%
  # quantile of all the 2.5% quantiles, and the 97.5% quantile of all the 
  # 97.5% quantiles.
  dplyr::summarize(
    dplyr::group_by(EpiEstim.result, date),
    Rt    = mean(I(mean)), # The 'I()' call preserves meaning of 'mean'
    Rt.lo = quantile(lo, 0.025),
    Rt.hi = quantile(hi, 0.975)
  ) -> final_values

  if (graph == FALSE) 
    return(final_values)

  ggplot2::ggplot(final_values, aes(x = date)) +
    geom_hline(
      yintercept = 1,
      color = "red",
      size = 0.5,
      show.legend = FALSE
    ) +
    geom_line(aes(y = Rt)) + 
    geom_ribbon(aes(y = Rt, ymin=Rt.lo, ymax=Rt.hi), alpha=0.3) +
    scale_x_date(date_breaks = '1 week',
                 date_labels = "%b %d",
                 minor_breaks = NULL) +
    scale_y_continuous(
      trans = "log10",
      limits = c(0, 8),
      breaks = c(0.5, 1, 2, 4, 8),
      minor_breaks = NULL,
      expand = c(0,0)
    ) +
    labs(
      x = NULL,
      y = "Rt", 
      title = "Effective Reproduction Number Estimate"
    ) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        size = rel(3/4), angle = 45, hjust = 1, vjust = 1
      )
    )
  # return(Rt_df)
}

#' @export
RtNaiveEstim <- function(cc,
                         window = 5, 
                         mean.si = 4.7,
                         std.si = 2.9,
                         graph = TRUE) {

  day_start <- 2
  day_end   <- day_start + window - 1
  n_days    <- cc$config$N_days
  ndb       <- cc$config$N_days_before

  first_date <- as.Date(cc$config$first_date, origin = '1970-01-01')
  cases <- cc$config$obs_cas

  EpiEstim::make_config(
    list(
      t_start = seq(day_start, n_days - window + 1), 
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
  result <- dplyr::transmute(result, day = t_start, y = `Median(R)`,
                      ymin = `Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`)

  if (graph == FALSE)
    return(
      dplyr::transmute(
        result,
        date = as.Date(first_date, origin = '1970-01-01') + lubridate::days(day - 1),
        NaiveRt = y, NaiveRt.lo = ymin, NaiveRt.hi = ymax
      )
    )

  p <- ggplot2::ggplot(
    result, aes(x = as.Date(first_date) + lubridate::days(day - 1))
  ) +
    geom_hline(
      yintercept = 1,
      color = "red",
      size = 0.5,
      show.legend = FALSE
    ) +
    geom_line(aes(y = y)) + 
    geom_ribbon(aes(y = y, ymin=ymin, ymax=ymax), alpha=0.3) +
    scale_x_date(date_breaks = '1 week',
                 date_labels = "%b %d",
                 minor_breaks = NULL,
                 limits = c(first_date, NA)) +
    coord_cartesian(ylim = c(0, 8)) +
    scale_y_continuous(
      trans = "log10",
      limits = c(0, 8),
      breaks = c(0.5, 1, 2, 4, 8),
      minor_breaks = NULL,
      expand = c(0,0)
    ) +
    labs(
      x = NULL,
      y = "Rt", 
      title = "Naive Effective Reproduction Number Estimate"
    ) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        size = rel(3/4), angle = 45, hjust = 1, vjust = 1
      )
    )

  p
}


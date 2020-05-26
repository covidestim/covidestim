#' @export
#' @rdname RtEst.covidcast_result
RtEst <- function(...) UseMethod('RtEst')

#' Calculate the Effective Reprouction Number from CovidCast output
#'
#' Creates a figure with Rt estimates over time and 95% CI and undelrying df. 
#'
#' @param cc The result of calling \code{\link{run}}. An object of class
#'   \code{covidcast_result}.
#'   
#' @param sample_fraction The fraction of iterations to sample for the Rt
#' calualation; defaults to 2/3 samples after warm-up.  
#' 
#' @param window The size of the moving window over which Rt should be 
#' estimated; defaults to 7 days. 
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
RtEst.covidcast_result <- function(cc, 
                                   sample_fraction = (2/3), 
                                   window = 7, 
                                   mean.si = 4.7,
                                   std.si = 2.9,
                                   pdf.name = "covidcast_Rt.pdf") {
  
  fit       <- cc$extracted 
  tot_iter  <- 500 # cc$iter
  warm      <- 400 # cc$warmup
  iter      <- (tot_iter - warm) * 3 # cc$chains
  n_sample  <- round(sample_fraction * iter)
  day_start <- 2
  day_end   <- window + day_start - 1
  n_days    <- cc$config$N_days
  ndb       <- cc$config$N_days_before

  sample_iter <- sample(iter, n_sample)

  inf <- as.data.frame(fit[["new_inf"]]) %>% 
    gather(key = day, value = I) %>%
    mutate(day = as.numeric(substr(day, start = 2, stop = 4))) %>%
    filter(day > ndb) %>%
    group_by(day) %>% 
    mutate(iter = 1:iter) %>% 
    ungroup() 

  # number of moving windows in which Rt will be estimated     
  # total day - window works if 1st day is day 2 of data
  windows <- length(unique(inf$day)) - window 

  # create matrices to store output of epi estim call 
  # where each row is value for a window, each column is estimates from an iter
  mn <- matrix(nrow = windows, ncol = n_sample)
  lo <- matrix(nrow = windows, ncol = n_sample)
  hi <- matrix(nrow = windows, ncol = n_sample)

  #' @import EpiEstim
  # now we generate estimates of the mean, upper, and lower bounds of Rt from 
  # each sampled iteration of case data. 
  for(i in seq_along(sample_iter)) {
    inf2 <- filter(inf, iter == sample_iter[i]) 
    
    EpiEstim::make_config(
      list(
        t_start = seq(day_start, nrow(inf2) - (day_end - day_start)), 
        t_end   = seq(day_end,   nrow(inf2)),
        mean_si = mean.si, 
        std_si  = std.si
      )
    ) -> config

    EpiEstim::estimate_R(
      inf2$I,
      method = "parametric_si",
      config = config
    ) -> out

    R_est <- out[["R"]]
    
    mn[,i] <- R_est$`Mean(R)`
    lo[,i] <- R_est$`Quantile.0.025(R)`
    hi[,i] <- R_est$`Quantile.0.975(R)`
  }

  # convert matriaces into df, calculate bounds on each estimate
  make_est_df <- function(x){
    as.data.frame(x) %>% 
    rowid_to_column() %>% 
    gather(key = "iter", 
           value = "est", 
           2:(windows+1)) %>%
    group_by(rowid) %>%
    summarise(mn = mean(est), 
              lo = quantile(est, 0.025), 
              hi = quantile(est, 0.975))
  }

  Rt    <- make_est_df(mn)
  lower <- make_est_df(lo) 
  upper <- make_est_df(hi) 

  Rt_df <- as.data.frame(cbind(Rt$mn, lower$lo, upper$hi)) %>% 
    mutate(day = seq(1,nrow(Rt),1)) 

  first_date <- as.Date(cc$config$first_date, origin = '1970-01-01')

  ggplot2::ggplot(
    Rt_df, aes(x = first_date + lubridate::days(day - 1))
  ) +
    geom_hline(
      yintercept = 1,
      color = "red",
      size = 0.5,
      show.legend = FALSE
    ) +
    geom_line(aes(y = V1)) + 
    geom_ribbon(aes(y = V1, ymin=V2, ymax=V3), alpha=0.3) +
    scale_x_date(date_breaks = '1 week',
                 date_labels = "%b %d",
                 minor_breaks = NULL) +
    scale_y_continuous(
      limits = c(0, 8),
      breaks = c(0, 1, 2, 4, 8),
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
RtNaiveEstim <- function(ccr) {

  first_date <- as.Date(ccr$config$first_date, origin = '1970-01-01')
  cases <- ccr$config$obs_cas

  EpiEstim::make_config(
    list(
      t_start = seq(8, length(cases) - 6),
      t_end   = seq(14, length(cases)),
      mean_si = 4.7,
      std_si  = 2.9
    )
  ) -> cfg

  EpiEstim::estimate_R(
    cases,
    method = "parametric_si",
    config = cfg
  ) -> Rt_Est

  result <- as_tibble(Rt_Est[["R"]])
  result <- transmute(result, day = t_start, y = `Median(R)`,
                      ymin = `Quantile.0.025(R)`, ymax=`Quantile.0.975(R)`)

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
      breaks = c(0, 1, 2, 4, 8),
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


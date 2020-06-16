#' @export
#' @rdname viz.covidestim_result
viz <- function(...) UseMethod('viz')

#' Visualize results of Covidestim
#'
#' Returns several graphs of Covidestim input and estimates
#'
#' @param ccr The result of calling \code{\link{run}}. An object of class
#'   \code{covidestim_result}.
#'
#' @param renderPDF If \code{FALSE}, figures are returned as a list of
#'   \code{ggplot} objects. If a string, a \code{pdf} will be rendered to that
#'   location, and and the list of \code{ggplot} objects will still be
#'   returned.
#'
#' @return A list with the following keys:
#'
#' \itemize{
#'   \item \code{observedVsFitted} A \code{ggplot} of input cases/deaths data, and the
#'     fit of the model to that data.
#'   \item \code{infectionsAndCases} A \code{ggplot} of estimated daily infections,
#'     compared to reported/observed cases.
#' }
#'
#' If \code{include.RtEstim == TRUE}:
#'
#' \itemize{
#'   \item \code{RtEstim} A \code{ggplot} object of estimated Rt
#'   \item \code{RtEstimNaive} A \code{ggplot} object of Rt estimated
#'     \strong{purely} from input data
#' }
#'
#' @seealso \code{\link{summary.covidestim_result}} for more details on the
#'   quantities being plotted
#'
#' @export
viz.covidestim_result <- function(ccr, include.RtEstim = TRUE, renderPDF = FALSE) {

  # Prep all the intermediate representations of the data that are ultimately
  # used to plot everything
  run_summary <- summary(ccr, include.RtEstim = include.RtEstim)

  first_date <- as.Date(ccr$config$first_date, origin = '1970-01-01')
  ndays      <- ccr$config$N_days


  # Kind of hack-y way to get the "input data frame"
  dplyr::bind_cols(
    date = seq(first_date, first_date + lubridate::days(ndays - 1), by = '1 day'),
    cases = ccr$config$obs_cas,
    deaths = ccr$config$obs_die
  ) -> input_data

  # Plot the two main graphs
  list(
    observedVsFitted = viz_observed_and_fitted(run_summary, input_data),
    infectionsAndCases = viz_all_cases_to_data(run_summary, input_data)
  ) -> result

  if (include.RtEstim)
    result <- append(result,
      list(
        RtEstim      = viz_RtEstim(run_summary),
        RtEstimNaive = viz_RtEstimNaive(run_summary)
      )
    )

  if (!identical(renderPDF, FALSE))
    ggsave(file = renderPDF, width = 8.5, height = 11) 

  result
}

viz_RtEstim <- function(run_summary) {

  ggplot2::ggplot(run_summary, aes(x = date)) +
    geom_hline(
      yintercept = 1,
      color = "red",
      size = 0.5,
      show.legend = FALSE
    ) +
    geom_line(aes(y = Rt), na.rm = TRUE) + 
    geom_ribbon(aes(y = Rt, ymin=Rt.lo, ymax=Rt.hi), alpha=0.3, na.rm = TRUE) +
    scale_x_date(date_breaks = '1 week',
                 date_labels = "%b %d",
                 minor_breaks = NULL) +
    scale_y_log10(breaks = c(seq(0.5, 1.5, 0.1), 1, 2, 3, 4, 5),
                  minor_breaks = NULL) +
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
}

viz_RtEstimNaive <- function(run_summary) {

  ggplot2::ggplot(run_summary, aes(x = date)) +
    geom_hline(
      yintercept = 1,
      color = "red",
      size = 0.5,
      show.legend = FALSE
    ) +
    geom_line(aes(y = NaiveRt), na.rm = TRUE) + 
    geom_ribbon(aes(y = NaiveRt, ymin=NaiveRt.lo, ymax=NaiveRt.hi), alpha=0.3,
                na.rm = TRUE) +
    scale_x_date(date_breaks = '1 week',
                 date_labels = "%b %d",
                 minor_breaks = NULL,
                 limits = c(min(run_summary$date), NA)) +
    scale_y_log10(breaks = c(seq(0.5, 1.5, 0.1), 1, 2, 3, 4, 5),
                  minor_breaks = NULL) +
    coord_cartesian(ylim = c(0, 8)) +
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
}

#' @import ggplot2
viz_observed_and_fitted <- function(run_summary, input_data) {

  first_date <- min(run_summary$date)

  custom_logscale <- scale_y_continuous(
    name = "Count (log scale)",
    trans = "log1p",
    labels = scales::label_number_si(),
    breaks = function(l) c(0, scales::breaks_log(n=5)(c(1, l[2]))),
    minor_breaks = NULL
  )

  ggplot2::ggplot(run_summary, aes(x = date)) + 
    geom_point(data = input_data, aes(y = cases, color = "ocas"),  size = rel(0.6)) + 
    geom_line(data  = input_data, aes(y = cases, color = "ocas")) + 

    geom_point(data = input_data, aes(y = deaths, color = "odth"), size = rel(0.6)) + 
    geom_line(data  = input_data, aes(y = deaths, color = "odth")) + 

    geom_line(aes(y = cases.fitted, color = "fcas")) + 
    geom_ribbon(
      aes(ymin = cases.fitted.lo, ymax = cases.fitted.hi, color = "fcas"),
      alpha = 0.3, linetype = 'dashed'
    ) +

    geom_line(aes(y = deaths.fitted)) + 
    geom_ribbon(
      aes(ymin = deaths.fitted.lo, ymax = deaths.fitted.hi, color = "fdth"),
      alpha = 0.3, linetype = 'dashed'
    ) +

    custom_logscale +
    scale_color_manual(
      values = c('#a6cee3','#b2df8a','#1f78b4','#33a02c'),
      labels = c("ocas" = "Observed Cases",
                 "odth" = "Observed Deaths",
                 "fcas" = "Fitted Cases",
                 "fdth" = "Fitted Deaths"),
      guide = guide_legend(
        override.aes = list(
          linetype = c(rep("solid", 2), rep("dashed", 2))
        ),
        nrow = 2
      )
    ) +
    scale_x_date(date_breaks = '1 week',
                 date_labels = "%b %d",
                 minor_breaks = NULL) +
    labs(x = NULL,
         title = "Observed and Fitted Cases and Deaths",
         color = "") +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        size = rel(3/4), angle = 45, hjust = 1, vjust = 1
      ),
      legend.position = "bottom",
      legend.text = element_text(size=rel(0.75))
    )
}

#' @import ggplot2
viz_all_cases_to_data <- function(run_summary, input_data) {

  ggplot2::ggplot(run_summary, aes(x = date)) + 
    geom_point(data = input_data,
               aes(y = cases, color = "ocas"),
               size = rel(0.6)) +
    geom_line(data = input_data,
              aes(y = cases, color = "ocas")) +
    geom_line(aes(y = infections, color = "ninf")) + 
    geom_ribbon(
      aes(ymin = infections.lo, ymax = infections.hi),
      alpha=0.3,
      linetype = 2
    ) +
    scale_color_manual(
      values = c('#a6cee3','#fc8d62'),
      labels = c("ninf" = "New Infections", "ocas" = "Observed Cases")
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      labels = scales::label_number_si(),
      breaks = c(0, 10^(1:6)),
      minor_breaks = NULL,
      limits = c(0, NA)
    ) +
    scale_x_date(date_breaks = '1 week',
                 date_labels = "%b %d",
                 minor_breaks = NULL) +
    labs(
      x = NULL,
      y = "Count",
      title = "Modeled New Infections",
      color = NULL
    ) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        size = rel(3/4), angle = 45, hjust = 1, vjust = 1
      ),
      legend.position = "bottom",
      legend.text = element_text(size=rel(0.75))
    )
}

# comparison: fitted, diagnosed, all cases

#' @import ggplot2
viz_modeled_cases <- function(cc, fit_to_data, diag, deltas) {

  first_date <- as.Date(cc$config$first_date, origin = '1970-01-01')

  ggplot2::ggplot(
    mapping = aes(x     = first_date + lubridate::days(day - 1),
                  y     = median,
                  color = outcome)
  ) + 
    geom_line(data = filter(fit_to_data, outcome == "occur_cas")) + 
    geom_ribbon(
      data = filter(fit_to_data, outcome == "occur_cas"),
      aes(ymin=lo, ymax=hi),
      linetype = 2,
      alpha=0.2
    ) +
    geom_line(data = filter(diag, outcome == "diag_all")) +
    geom_ribbon(
      data = filter(diag, outcome == "diag_all"),
      aes(ymin=lo, ymax=hi),
      linetype = 2,
      alpha=0.2
    ) +
    # geom_line(data = filter(deltas, outcome == "new_inf")) + 
    # geom_ribbon(
    #   data = filter(deltas, outcome == "new_inf"),
    #   aes(ymin=lo, ymax=hi),
    #   linetype = 2,
    #   alpha=0.2
    # ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      labels = scales::label_number_si(),
      # breaks = function(lims) c(0, scales::breaks_log(n = 5)(1, lims[2])),
      breaks = 10^(1:6),
      minor_breaks = NULL,
      limits = c(0, NA)
    ) +
    scale_x_date(date_breaks = '1 week',
                 date_labels = "%b %d",
                 minor_breaks = NULL) +
    scale_color_manual(
      values = c('#fc8d62','#66c2a5','#8da0cb'), 
      breaks = c("diag_all", "occur_cas"),
      labels = c("diag_all" = "Modeled Diagnosed Cases", 
                 "occur_cas" = "Fitted Reported Cases")) +
    # scale_color_manual(
    #   values = c('#fc8d62','#66c2a5','#8da0cb'), 
    #   breaks = c("new_inf", "diag_all", "occur_cas"),
    #   labels = c("new_inf" = "Modeled New Infections",
    #              "diag_all" = "Modeled Diagnosed Cases", 
    #              "occur_cas" = "Fitted Reported Cases")) +
    labs(x = NULL,
         y = "Count",
         # title = "Modeled New Infections, Diagnosed Cases, and Reported COVID-19 Cases",
         title = "Modeled Diagnoses and Case Reports",
         color = NULL) +
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        size = rel(3/4), angle = 45, hjust = 1, vjust = 1
      ),
      legend.position = "bottom",
      legend.text = element_text(size=rel(0.75))
    )
}


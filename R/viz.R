#' @export
viz <- function(...) UseMethod('viz')

#' Visualize results of a Covidestim run
#'
#' Returns two graphs of Covidestim input and estimates
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
#' @seealso \code{\link{summary.covidestim_result}} for more details on the
#'   quantities being plotted
#'
#' @export
viz.covidestim_result <- function(ccr, renderPDF = FALSE) {

  # Prep all the intermediate representations of the data that are ultimately
  # used to plot everything
  run_summary <- summary(ccr)

  first_date <- as.Date(ccr$config$first_date, origin = '1970-01-01')
  nweeks      <- ccr$config$N_weeks


  # Kind of hack-y way to get the "input data frame"
  dplyr::bind_cols(
    date = seq(first_date, first_date + lubridate::weeks(nweeks - 1), by = '1 week'),
    cases = ccr$config$obs_cas,
    hosp = ccr$config$obs_hosp,
    deaths = ccr$config$obs_die,
  ) -> input_data

  # Plot the two main graphs
  list(
    observedVsFitted = viz_observed_and_fitted(run_summary, input_data),
    infectionsAndCases = viz_all_cases_to_data(run_summary, input_data)
  ) -> result

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
    geom_line(aes(y = r_t), na.rm = TRUE) + 
    geom_ribbon(aes(y = r_t, ymin=r_t_p2_5, ymax=r_t_p97_5), alpha=0.3, na.rm = TRUE) +
    geom_ribbon(aes(y = r_t, ymin=r_t_p25, ymax=r_t_p75), alpha=0.3, na.rm = TRUE) +
    scale_x_date(date_breaks = '1 month',
                 date_labels = "%b",
                 minor_breaks = NULL) +
    scale_y_log10(breaks = c(0.5, 0.7, 1, 1.5, 2, 3, 4, 5),
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

#' @import ggplot2
viz_observed_and_fitted <- function(run_summary, input_data) {

  first_date <- min(run_summary$date)

  custom_logscale <- scale_y_continuous(
    name = "Count (log scale)",
    trans = "log1p",
    labels = scales::label_number(scale_cut = scales::cut_short_scale()),
    breaks = function(l) c(0, scales::breaks_log(n=5)(c(1, l[2]))),
    minor_breaks = NULL
  )

  ggplot2::ggplot(run_summary, aes(x = date)) + 
    geom_point(data = input_data, aes(y = cases, color = "ocas"),  size = rel(0.6)) + 
    geom_line(data  = input_data, aes(y = cases, color = "ocas")) + 

    geom_point(data = input_data, aes(y = hosp, color = "ohsp"), size = rel(0.6)) + 
    geom_line(data  = input_data, aes(y = hosp, color = "ohsp")) + 

    geom_point(data = input_data, aes(y = deaths, color = "odie"), size = rel(0.6)) + 
    geom_line(data  = input_data, aes(y = deaths, color = "odie")) + 

    geom_line(aes(y = fitted_cases, color = "fcas")) + 
    geom_ribbon(
      aes(ymin = fitted_cases_p2_5, ymax = fitted_cases_p97_5, color = "fcas"),
      alpha = 0.3, linetype = 'dashed'
    ) +
    geom_ribbon(
      aes(ymin = fitted_cases_p25, ymax = fitted_cases_p75, color = "fcas"),
      alpha = 0.3, linetype = 'dashed'
    ) +

    geom_line(aes(y = fitted_hospitalizations)) + 
    geom_ribbon(
      aes(ymin =fitted_hospitalizations_p2_5, ymax = fitted_hospitalizations_p97_5, color = "fhsp"),
      alpha = 0.3, linetype = 'dashed'
    ) +
    geom_ribbon(
      aes(ymin =fitted_hospitalizations_p25, ymax = fitted_hospitalizations_p75, color = "fhsp"),
      alpha = 0.3, linetype = 'dashed'
    ) +

    geom_line(aes(y = fitted_deaths)) + 
    geom_ribbon(
      aes(ymin =fitted_deaths_p2_5, ymax = fitted_deaths_p97_5, color = "fdie"),
      alpha = 0.3, linetype = 'dashed'
    ) +
    geom_ribbon(
      aes(ymin =fitted_deaths_p25, ymax = fitted_deaths_p75, color = "fdie"),
      alpha = 0.3, linetype = 'dashed'
    ) +
    
    custom_logscale +
    scale_color_manual(
      values = c('#a6cee3','#b2df8a','#fdbb84','#1f78b4','#33a02c', '#e34a33'),
      labels = c("ocas" = "Observed Cases",
                 "ohsp" = "Observed Hospitalizations",
                 "odie" = "Observed Deaths",
                 "fcas" = "Fitted Cases",
                 "fhsp" = "Fitted Hospitalizations",
                 "fdie" = "Fitted Deaths"),
      guide = guide_legend(
        override.aes = list(
          linetype = c(rep("solid", 3), rep("dashed", 3))
        ),
        nrow = 2
      )
    ) +
    scale_x_date(date_breaks = '1 month',
                 date_labels = "%b",
                 minor_breaks = NULL) +
    labs(x = NULL,
         title = "Observed and Fitted Cases, Deaths and Hospitalizations",
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
      aes(ymin = infections_p2_5, ymax = infections_p97_5),
      alpha=0.3,
      linetype = 2
    ) +
    geom_ribbon(
      aes(ymin = infections_p25, ymax = infections_p75),
      alpha=0.3,
      linetype = 2
    ) +
    scale_color_manual(
      values = c('#a6cee3','#fc8d62'),
      labels = c("ninf" = "New Infections", "ocas" = "Observed Cases")
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      breaks = c(0, 10^(1:6)),
      minor_breaks = NULL,
      limits = c(0, NA)
    ) +
    scale_x_date(date_breaks = '1 month',
                 date_labels = "%b",
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
  nweeks      <- ccr$config$N_weeks
  
  ggplot2::ggplot(
    mapping = aes(x     = first_date + lubridate::days(day - 1),
                  y     = median,
                  color = outcome)
  ) + 
    geom_line(data = filter(fit_to_data, outcome == "fitted_cases")) + 
    geom_ribbon(
      data = filter(fit_to_data, outcome == "fitted_cases"),
      aes(ymin=lo, ymax=hi),
      linetype = 2,
      alpha=0.2
    ) +
    geom_line(data = filter(diag, outcome == "diagnoses")) +
    geom_ribbon(
      data = filter(diag, outcome == "diagnoses"),
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
      labels = scales::label_number(scale_cut = scales::cut_short_scale()),
      # breaks = function(lims) c(0, scales::breaks_log(n = 5)(1, lims[2])),
      breaks = 10^(1:6),
      minor_breaks = NULL,
      limits = c(0, NA)
    ) +
    scale_x_date(date_breaks = '1 month',
                 date_labels = "%b",
                 minor_breaks = NULL) +
    scale_color_manual(
      values = c('#fc8d62','#66c2a5','#8da0cb'), 
      breaks = c("diagnoses", "fitted_cases"),
      labels = c("diagnoses" = "Modeled Diagnosed Cases", 
                 "fitted_cases" = "Fitted Reported Cases")) +
    # scale_color_manual(
    #   values = c('#fc8d62','#66c2a5','#8da0cb'), 
    #   breaks = c("new_inf", "diagnoses", "fitted_cases"),
    #   labels = c("new_inf" = "Modeled New Infections",
    #              "diagnoses" = "Modeled Diagnosed Cases", 
    #              "fitted_cases" = "Fitted Reported Cases")) +
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


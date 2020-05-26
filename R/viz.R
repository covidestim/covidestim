#' @export
#' @rdname viz.covidcast_result
viz <- function(...) UseMethod('viz')

#' Visualize results of Covidcast
#'
#' Creates several graphs of covidcast results
#'
#' @param cc The result of calling \code{\link{run}}. An object of class
#'   \code{covidcast_result}.
#'
#' @param renderPDF A logical scalar. If true, will render a pdf to the working
#'   directory with name \code{covidcast_figures.pdf}.
#'
#' @return At this time, returns a \code{\link[ggplot]{ggplot}} object with the
#'   last plot rendered. This is likely to change.
#'
#' @export
viz.covidcast_result <- function(cc, renderPDF = FALSE) {

  pdfname <- glue("covidcast-output.pdf")

  # if (renderPDF)
  #   pdf(file = pdfname, width = 8, height = 6) 

  # Prep all the intermediate representations of the data that are ultimately
  # used to plot everything
  intermediate_objects <- viz_prep(cc)

  input_data  <- intermediate_objects$input_data
  fit_to_data <- intermediate_objects$fit_to_data 
  diag        <- intermediate_objects$diag 
  deltas      <- intermediate_objects$deltas 

  N_days_before <- cc$config$N_days_before

  # if (renderPDF) {
  #   dev.off()
  #   message(glue("PDF successfully rendered to ./{pdfname}"))
  # }

  # Plot the first four graphs. The remaining graphs aren't active yet.
  list(
    p1 = viz_observed_and_fitted(cc, input_data, fit_to_data, N_days_before),
    p2 = viz_all_cases_to_data(cc, input_data, deltas),
    p3 = viz_modeled_cases(cc, fit_to_data, diag, deltas)
#     p4 = viz_incidence(cc, fit_to_data, diag, deltas)
  ) 
}


#' COVID viz 
#' @importFrom dplyr mutate group_by summarise ungroup arrange
#' @importFrom magrittr %>%
viz_prep <- function(obj) {
  result  <- obj$result
  fit     <- obj$extracted
  datList <- obj$config

  summary_fit <- rstan::summary(result)

  N_days_before <- obj$config$N_days_before

  outcomeGen <- function(idx) {
    as.data.frame(fit[[idx]]) %>% 
    tidyr::gather(key = day, value = estim) %>%
    mutate(outcome = idx)
  }

  new_inf   <- outcomeGen("new_inf")
  new_sym   <- outcomeGen("new_sym")
  new_sev   <- outcomeGen("new_sev")
  new_die   <- outcomeGen("new_die")

  new_sym_dx  <- outcomeGen("new_sym_dx")
  new_sev_dx  <- outcomeGen("new_sev_dx")
  new_die_dx  <- outcomeGen("new_die_dx")
  diag_all    <- outcomeGen("diag_all")

  occur_cas <- outcomeGen("occur_cas")
  occur_die <- outcomeGen("occur_die")

  reformat_staninputs <- function(vec, outcome)
    tibble::as_tibble(list(estim = vec)) %>%
      mutate(day = 1:n(), outcome = outcome)

  obs_cas <- reformat_staninputs(obj$config$obs_cas, "obs_cas")
  obs_die <- reformat_staninputs(obj$config$obs_die, "obs_die")

  deltas <- rbind(new_inf, new_sym, new_sev, new_die) %>%
    group_by(day, outcome) %>%
      summarise(median = median(estim), 
      lo = quantile(estim, 0.025), 
      hi = quantile(estim, 0.975)) %>%
    ungroup() %>%
    mutate(day = as.numeric(substr(day, start = 2, stop = 4)) - N_days_before) %>%
    arrange(day)

  diag <- rbind(new_sym_dx, new_sev_dx, diag_all, new_die_dx) %>%
    group_by(day, outcome) %>%
    summarise(median = median(estim), 
              lo = quantile(estim, 0.025), 
              hi = quantile(estim, 0.975)) %>%
    ungroup() %>%
    mutate(day = as.numeric(substr(day, start = 2, stop = 4)) - N_days_before) %>%
    arrange(day)


  fit_to_data <- rbind(occur_cas, occur_die) %>%
    group_by(day, outcome) %>%
    summarise(median = median(estim), 
              lo = quantile(estim, 0.025), 
              hi = quantile(estim, 0.975)) %>%
    ungroup() %>%
    mutate(day = as.numeric(substr(day, start = 2, stop = 4)) - N_days_before) %>%
    arrange(day)


  input_data <- rbind(obs_cas, obs_die) %>%
                group_by(day, outcome) %>% summarise(median = median(estim), 
                                         lo = quantile(estim, 0.025), 
                                         hi = quantile(estim, 0.975)) %>%
                ungroup() %>% arrange(day)

  list(deltas=deltas, diag=diag, fit_to_data=fit_to_data, input_data=input_data)
}


#' @import ggplot2
viz_observed_and_fitted <- function(cvc_result, input_data, fit_to_data,
                                    N_days_before) {

  first_date <- as.Date(cvc_result$config$first_date, origin = '1970-01-01')

  over1k_old <- scale_y_continuous(
    name = "Count (log1p scale)",
    trans = scales::pseudo_log_trans(base = 10),
    labels = scales::label_number_si(),
    breaks = function(lims) sort(c(0, scales::breaks_log(n = 5)(0.01, lims[2]))),
    minor_breaks = NULL
  )

  over1k <- scale_y_continuous(
    name = "Count (log scale)",
    trans = "log1p",
    labels = scales::label_number_si(),
    breaks = function(l) c(0, scales::breaks_log(n=5)(c(1, l[2]))),
    minor_breaks = NULL
  )
  
  under1k <- scale_y_continuous(name = "Count")

  # dynamicScale <- if (max(input_data$median) > 1e3) over1k else under1k

  ggplot2::ggplot(input_data,
                  aes(x = first_date + lubridate::days(day - 1),
                      y = median,
                      color = outcome)) + 
    geom_point(data = input_data, size = rel(0.6)) + 
    geom_line(data  = input_data) + 
    geom_line(data  = filter(fit_to_data, day > 0)) + 
    geom_ribbon(
      data = filter(fit_to_data, day > 0),
      aes(ymin = lo, ymax = hi),
      alpha = 0.3,
      linetype = 'dashed'
    ) +
    # dynamicScale +
    over1k +
    scale_color_manual(
      values = c('#a6cee3','#b2df8a','#1f78b4','#33a02c'),
      labels = c("Observed Cases",
                 "Observed Deaths",
                 "Fitted Cases",
                 "Fitted Deaths"),
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

  ## all cases to data 

#' @import ggplot2
viz_all_cases_to_data <- function(cc, input_data, deltas) {

  first_date <- as.Date(cc$config$first_date, origin = '1970-01-01')

  ggplot2::ggplot(mapping = aes(x = first_date + lubridate::days(day - 1),
                                y = median,
                                color = outcome)) + 
    geom_point(data = filter(input_data, outcome == "obs_cas"),
               size = rel(0.6)) +
    geom_line(
      data = filter(input_data, outcome == "obs_cas")
    ) + 
    geom_line(data = filter(deltas, outcome == "new_inf", day >= 1)) + 
    geom_ribbon(
      data = filter(deltas, outcome == "new_inf", day >= 1),
      aes(ymin=lo, ymax=hi),
      alpha=0.3,
      linetype = 2
    ) +
    scale_color_manual(
      values = c('#a6cee3','#fc8d62'),
      labels = c("New Infections", "Observed Cases")
    ) +
    scale_y_continuous(
      trans = scales::pseudo_log_trans(base = 10),
      labels = scales::label_number_si(),
      # breaks = function(lims) c(0, scales::breaks_log(n = 5)(1, lims[2])),
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

# plotting: all, diagnosed, reported outcomes

#' @import ggplot2
viz_incidence <- function(cc, fit_to_data, diag, deltas) {
  ggplot2::ggplot(
    mapping = aes(x = day, y = median, color = outcome)
  ) +
    geom_line(data = fit_to_data) +
    geom_line(data = diag) +
    geom_line(data = deltas) +
    scale_y_continuous(trans = "log10") +
    scale_color_manual(
      values = c('#08519c','#3182bd', '#bdd7e7', '#238b45',
                 '#74c476','#bae4b3', '#6a51a3', '#cbc9e2'
                 # 'green', 'purple', 'orange', 'blue'
                 ), 
      breaks = c("new_inf",
                 "new_sym",
                 #"new_sev",
                 "new_die",
                 "diag_all",
                 "new_sym_dx",
                 #"new_sev_dx",
                 "new_die_dx",
                 "occur_cas",
                 "occur_die",
                 #"obs_cas", 
                 #"obs_die"
                 ),
      labels = c("Modeled New Infections",
                 "Modeled Symptomatic Cases", 
                 #"Modeled New Severe Cases",
                 "Modeled Deaths",
                 "Modeled All Diagnosed", 
                 #"Modeled Diagnosed at Symptomatic", 
                 "Modeled Diagnosed at Severe", 
                 "Modeled Diagnosed at Death", 
                 "Fitted Cases",
                 "Fitted Deaths"
                 #"Observed Cases",
                 #"Observed Deaths"
                 )
    ) +
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
    labs(x = NULL,
         y = "Count",
         title = "Median Modeled Outcomes Compared to Data",
         color = "")
    theme_linedraw() +
    theme(
      axis.text.x = element_text(
        size = rel(3/4), angle = 45, hjust = 1, vjust = 1
      ),
      legend.position = "bottom",
      legend.text = element_text(size = rel(0.75))
    )
}

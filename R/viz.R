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

  if (renderPDF)
    pdf(file = pdfname, width = 8, height = 6) 

  # Prep all the intermediate representations of the data that are ultimately
  # used to plot everything
  intermediate_objects <- viz_prep(cc)

  input_data  <- intermediate_objects$input_data
  fit_to_data <- intermediate_objects$fit_to_data 
  diag        <- intermediate_objects$diag 
  deltas      <- intermediate_objects$deltas 

  N_days_before <- cc$config$N_days_before

  # Plot the first four graphs. The remaining graphs aren't active yet.
  print(viz_observed_and_fitted(input_data, fit_to_data, N_days_before))
  print(viz_all_cases_to_data(input_data, deltas))
  print(viz_modeled_cases(fit_to_data, diag, deltas))
  print(viz_incidence(fit_to_data, diag, deltas))

  if (renderPDF) {
    dev.off()
    message(glue("PDF successfully rendered to ./{pdfname}"))
  }
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


  ## PLOTS
  ## fitted to data 
  # pdfnam <- "Figures_covidcast_4.21.pdf"
  # pdf(file=pdfnam,width=8, height=6) 


  # needs: input_data, fit_to_data

#' @import ggplot2
viz_observed_and_fitted <- function(input_data, fit_to_data,
                                    N_days_before) {
  ggplot2::ggplot() + 
    geom_point(
      data = input_data,
      aes(x = day, y = median, color = outcome)
    ) + 
    geom_line(
      data = input_data,
      aes(x = day, y = median, color = outcome),
      linetype = 2
    ) + 
    geom_line(
      data = filter(fit_to_data, day >=0),
      aes(x = day, y = median, color = outcome)
    ) + 
    geom_ribbon(
      data = filter(fit_to_data, day >= 0), 
      aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome),
      alpha=0.3
    ) +
    scale_color_manual(
      values = c('#a6cee3','#b2df8a','#1f78b4','#33a02c'),
      labels = c("Reported Cases", "Reported Deaths",
                "Fitted Reported Cases","Fitted Reported Deaths"),
      guide = guide_legend(override.aes = list(linetype = c(rep("solid", 2), rep("dashed", 2))))
    ) +
    labs(x = "Days Since March 1, 2020",
         y = "Count",
         title = "Observed and Fitted COVID-19 Cases, Deaths in NYC: March 2 - April 27",
         caption = "Data accessed April 28, 2020",
         color = "")
}

  ## all cases to data 

#' @import ggplot2
viz_all_cases_to_data <- function(input_data, deltas) {
  ggplot2::ggplot() + 
    geom_point(
      data = filter(input_data, outcome == "obs_cas"),
      aes(x = day, y = median, color = outcome)
    ) + 
    geom_line(
      data = filter(input_data, outcome == "obs_cas"),
      aes(x = day, y = median, color = outcome),
      linetype = 2
    ) + 
    geom_line(
      data = filter(deltas, outcome == "new_inf", day >=0),
      aes(x = day, y = median, color = outcome)
    ) + 
    geom_ribbon(
      data = filter(deltas, outcome == "new_inf", day >= 0),
      aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome),
      alpha=0.3
    ) +
    scale_color_manual(
      values = c('#a6cee3','#fc8d62'),
      labels = c("Modeled Cases", "Observed Cases"),
      guide = guide_legend(override.aes = list(linetype = c("dashed","solid")))
    ) +
    scale_y_log10() +
    labs(
      x = "Days Since Start",
      y = "Count",
      title = "Observed and Fitted COVID-19 Cases, Hosptialization, Deaths: NYC March 2 - April 27",
      subtitle = "with 95% uncertainty intervals",
      caption = "Data accessed April 28, 2020",
      color = ""
    )
}

# comparison: fitted, diagnosed, all cases

#' @import ggplot2
viz_modeled_cases <- function(fit_to_data, diag, deltas) {
  ggplot2::ggplot() + 
    geom_line(
      data = filter(fit_to_data, outcome == "occur_cas"),
      aes(x = day, y = median, color = outcome)) + 
    geom_ribbon(
      data = filter(fit_to_data, outcome == "occur_cas"),
      aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome),
      linetype = 2,
      alpha=0.2
    ) +
    geom_line(
      data = filter(diag, outcome == "diag_all"),
      aes(x = day, y = median, color = outcome)
    ) +
    geom_ribbon(
      data = filter(diag, outcome == "diag_all"),
      aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome),
      linetype = 2,
      alpha=0.2
    ) +
    geom_line(
      data = filter(deltas, outcome == "new_inf"),
      aes(x = day, y = median, color = outcome)) +
    geom_ribbon(
      data = filter(deltas, outcome == "new_inf"),
      aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome),
      linetype = 2,
      alpha=0.2
    ) +
    scale_y_log10() +
    scale_color_manual(
      values = c('#fc8d62','#66c2a5','#8da0cb'), 
      breaks = c("new_inf", "diag_all", "occur_cas"),
      labels = c("new_inf" = "Modeled New Infections",
                 "diag_all" = "Modeled Diagnosed Cases", 
                 "occur_cas" = "Fitted Reported Cases")) +
    labs(x = "Days Since March 1, 2020",
         y = "Count",
         title = "Modeled Infections, Diagnosed, and Reported COVID-19 Cases: NYC March 2 - April 27",
         subtitle = "with 95% uncertainty intervals",
         caption = "Data accessed April 28, 2020",
         color = "")
}

# plotting: all, diagnosed, reported outcomes

#' @import ggplot2
viz_incidence <- function(fit_to_data, diag, deltas) {
  ggplot2::ggplot() +
    geom_line(data = fit_to_data, 
              aes(x = day, y = median, color = outcome)) + 
    geom_line(data = diag, 
              aes(x = day, y = median, color = outcome)) +
    geom_line(data = deltas, 
              aes(x = day, y = median, color = outcome)) + 
    scale_y_log10() +
    scale_color_manual(
      values = c('#08519c','#3182bd', '#bdd7e7', '#238b45',
                 '#74c476','#bae4b3', '#6a51a3', '#cbc9e2',
                 'green', 'purple', 'orange', 'blue'), 
      breaks = c("new_inf",
                 "new_sym",
                 "new_sev",
                 "new_die",
                 "diag_all",
                 "new_sym_dx",
                 "new_sev_dx",
                 "new_die_dx",
                 "occur_cas",
                 "occur_die",
                 "obs_cas",
                 "obs_die"),
      labels = c("Modeled New Infections",
                 "Modeled Symptomatic Cases", 
                 "Modeled New Severe Cases",
                 "Modeled Deaths",
                 "Modeled All Diagnosed", 
                 "Modeled Diagnosed at Symptomatic", 
                 "Modeled Diagnosed at Severe", 
                 "Modeled Diagnosed at Death", 
                 "Fitted Cases",
                 "Fitted Deaths",
                 "Observed Cases",
                 "Observed Deaths")) +
    labs(x = "Days Since March 1, 2020",
         y = "Count",
         title = "Incidence Outcomes Compared to Data: NYC March 2 - April 27",
         subtitle = "with 95% uncertainty intervals",
         caption = "Data accessed April 28, 2020",
         color = "")
}


#  new.par <- par(mfrow=c(2,3))
#  mod <- fit[["p_sym_if_inf"]]
#  mn <- as.numeric(datList[["pri_p_sym_if_inf_a"]])
#  ss <- as.numeric(datList[["pri_p_sym_if_inf_b"]])
#  hist(mod, main = "P(Progression to Symptomatic)", xlab = "", ylab = "")
#  zz <- range(c(mod,0,1))
#  x_pts <- seq(zz[1],zz[2],length.out=101)
#  y_pts <- dbeta(x_pts,mn,ss)
#  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
#  
#  mod <- fit[["p_hos_if_sym"]]
#  mn <- as.numeric(datList[["pri_p_hos_if_sym_a"]])
#  ss <- as.numeric(datList[["pri_p_hos_if_sym_b"]])
#  hist(mod, main = "P(Progression to Hospitalizable)", xlab = "", ylab = "")
#  zz <- range(c(mod,0,1))
#  x_pts <- seq(zz[1],zz[2],length.out=101)
#  y_pts <- dbeta(x_pts,mn,ss)
#  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
#  
#  mod <- fit[["p_die_if_hos"]]
#  mn <- as.numeric(datList[["pri_p_die_if_hos_a"]])
#  ss <- as.numeric(datList[["pri_p_die_if_hos_b"]])
#  hist(mod, main = "P(Progression to Death)", xlab = "", ylab = "")
#  zz <- range(c(mod,0,1))
#  x_pts <- seq(zz[1],zz[2],length.out=101)
#  y_pts <- dbeta(x_pts,mn,ss)
#  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
#  
#  mod <- fit[["p_diag_if_inf"]]
#  mn <- as.numeric(datList[["pri_p_diag_if_inf_a"]])
#  ss <- as.numeric(datList[["pri_p_diag_if_inf_b"]])
#  hist(mod, main = "P(Diagnosis at Infectious)", xlab = "", ylab = "")
#  zz <- range(c(mod,0,1))
#  x_pts <- seq(zz[1],zz[2],length.out=101)
#  y_pts <- dbeta(x_pts,mn,ss)
#  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
#  
#  mod <- fit[["p_diag_if_sym"]]
#  mn <- as.numeric(datList[["pri_p_diag_if_sym_a"]])
#  ss <- as.numeric(datList[["pri_p_diag_if_sym_b"]])
#  hist(mod, main = "P(Diagnosis at Symptomatic)", xlab = "", ylab = "")
#  zz <- range(c(mod,0,1))
#  x_pts <- seq(zz[1],zz[2],length.out=101)
#  y_pts <- dbeta(x_pts,mn,ss)
#  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
#  
#  mod <- fit[["p_diag_if_hos"]]
#  mn <- as.numeric(datList[["pri_p_diag_if_hos_a"]])
#  ss <- as.numeric(datList[["pri_p_diag_if_hos_b"]])
#  hist(mod, main = "P(Diagnosis at Hospitalizable)", xlab = "", ylab = "")
#  zz <- range(c(mod,0,1))
#  x_pts <- seq(zz[1],zz[2],length.out=101)
#  y_pts <- dbeta(x_pts,mn,ss)
#  lines(x_pts,y_pts,col=mTrsp(4,200),lwd=2)
#
#  dev.off();

#
#  library("EpiEstim")
#  library("lubridate")
#
#  data1 <- group_by(new_inf, day) %>% summarise(cases = median(estim)) %>% 
#    mutate(day = substr(day, start = 2, stop =4)) %>%
#    arrange(as.numeric(day)) %>%
#    mutate(I = round(cases), 
#           dates = seq(as.Date("2020/02/21"), as.Date("2020/04/20"), by =1))
#
#  # assuming a mean SI of 4.1 days with a sd of 2.1 days
#  lognorm <- c(rlnorm(n = 100000, meanlog = 1.41, sdlog = 0.75))
#  lognorm2 <- table(round(lognorm))
#  lognorm_dist <- c(0, prop.table(lognorm2))
#
#  config <- make_config(list(t_start = seq(8, nrow(data1)-4), 
#                             t_end = seq(11, nrow(data1)-1),
#                             si_distr = lognorm_dist))
#
#  R_est  <- estimate_R(data1,
#                       method = "non_parametric_si",
#                       config = config)
#
#  print(estimate_R_plots(R_est, what = "R") + 
#    labs(title = "Preliminary R Effective Estimate, NYC", 
#         subtitle = "Assuming mean SI of 4.1 days (sd 2.1 days)") )
#
#  dev.off()
#  ############
#
#
#  ## comparison: diagnosed to fitted 
#  ggplot2::ggplot() + 
#    geom_line(data = diag, aes(x = day, y = median, color = outcome)) + 
#    geom_line(data = filter(fit_to_data, outcome != "occur_die"), 
#              aes(x = day, y = median, color = outcome), linetype = 2) + 
#    scale_color_manual(values = c('#a6cee3','#fb9a99', '#cab2d6', '#fdbf6f', '#1f78b4','#e31a1c'), 
#                       labels = c("Modeled All Diagnosed", "Modeled Diagnosed at Hospitalization", 
#                                  "Modeled Diagnoed at Infectious", "Modeled Diagnosed at Symptomatic",
#                                  "Fitted Cases", "Fitted Hospitalizations"), 
#                       guide = guide_legend(override.aes = list(
#                         linetype = c(rep("solid", 4), rep("dashed", 2))))) +
#    labs(x = "Days Since March 2, 2020",
#         y = "Count",
#         title = "Modeled (Diagnosed) and Fitted COVID-19 Cases in NYC, March 2 - April 15",
#         color = "")
#
#  ## comparison: modeled (all) to fitted input 
#  ggplot2::ggplot() + 
#    geom_line(data = filter(new, outcome != "new_sym"), aes(x = day, y = median, color = outcome)) + 
#    geom_line(data = fit_to_data, aes(x = day, y = median, color = outcome), linetype = 2) + 
#    scale_color_manual(values = c('#a6cee3','#b2df8a','#fb9a99','#1f78b4','#33a02c','#e31a1c'), 
#                       labels = c("Modeled All Deaths", "Modeled All Hospitalizations", "Modeled All Infections", 
#                                  "Fitted Cases", "Fitted Deaths", "Fitted Hospitalizations"), 
#                       guide = guide_legend(override.aes = list(
#                         linetype = c(rep("solid", 3), rep("dashed", 3))))) +
#    labs(x = "Days Since March 2, 2020",
#         y = "Count",
#         title = "Modeled (All) and Fitted COVID-19 Cases in NYC, March 2 - April 15",
#         color = "")
#
#  ###
#  ## SAME but with uncertainty 
#
#  ## comparison: diagnosed to fitted 
#  ggplot2::ggplot() + 
#    geom_line(data = diag, aes(x = day, y = median, color = outcome)) + 
#      geom_ribbon(data = diag,
#                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), 
#                linetype = 2, alpha=0.3) +
#    geom_line(data = filter(fit_to_data, outcome != "occur_die"), 
#              aes(x = day, y = median, color = outcome), linetype = 2) + 
#      geom_ribbon(data = filter(fit_to_data, outcome != "occur_die"), 
#                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), linetype = 2, alpha=0.3) +
#    scale_color_manual(values = c('#a6cee3','#fb9a99', '#cab2d6', '#fdbf6f', '#1f78b4','#e31a1c'), 
#                       labels = c("Modeled All Diagnosed", "Modeled Diagnosed at Hospitalization", 
#                                  "Modeled Diagnoed at Infectious", "Modeled Diagnosed at Symptomatic",
#                                  "Fitted Cases", "Fitted Hospitalizations"), 
#                       guide = guide_legend(override.aes = list(
#                         linetype = c(rep("solid", 4), rep("dashed", 2))))) +
#    labs(x = "Days Since March 2, 2020",
#         y = "Count",
#         title = "Modeled (Diagnosed) and Fitted COVID-19 Cases in NYC, March 2 - April 15",
#         color = "")
#
#  ## comparison: modeled (all) to fitted input 
#  ggplot2::ggplot() + 
#    geom_line(data = filter(new, outcome != "new_sym"), aes(x = day, y = median, color = outcome)) + 
#        geom_ribbon(data = filter(new, outcome != "new_sym"),
#                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), linetype = 2, alpha=0.3) +
#    geom_line(data = fit_to_data, aes(x = day, y = median, color = outcome)) + 
#        geom_ribbon(data = fit_to_data, 
#                aes(x = day, y = median, ymin=lo, ymax=hi, color = outcome), linetype = 2, alpha=0.3) +
#    scale_color_manual(values = c('#a6cee3','#b2df8a','#fb9a99','#1f78b4','#33a02c','#e31a1c'), 
#                       labels = c("Modeled All Deaths", "Modeled All Hospitalizations", "Modeled All Infections", 
#                                  "Fitted Cases", "Fitted Deaths", "Fitted Hospitalizations"), 
#                       guide = guide_legend(override.aes = list(
#                         linetype = c(rep("solid", 6))))) +
#    labs(x = "Days Since March 2, 2020",
#         y = "Count",
#         title = "Modeled (All) and Fitted COVID-19 Cases in NYC, March 2 - April 15",
#         color = "")

#' Comparison to data figure #1
#'
#' @param samps The result of \code{\link[rstan]{extract}}
#' @param datList The data block passed to \code{\link[rstan]{stan}}
#' @param diagData The fake diagnosis data returned from \code{\link{genFakeData}}
#'
#' @return Side effects of plotting
#' @export
viz_comparison_to_data_1 <- function(samps, datList, diagData) {
  ########### Comparison to data ###########
  # pdfnam <- "covidcast_comparison_to_data_1.pdf"
  # pdf(file = pdfnam, width = 9, height = 6.5)

  # N_days <- min(datList[["N_days"]], dim(diagData$rep_tri_conf_cases)[1])
  N_days <- datList[["N_days"]]
  print(N_days)
  N_days_tot <- datList[["N_days"]] + datList[["N_days_extra"]]

  max_delay <- max(diagData$days_delay) + 1
  cls <- colorRampPalette(c("red4", "red", "blue", "blue4"))(max_delay)

  otcm_mu <- apply(samps[["rep_tri_conf_cases_mu"]], 2:3, mean)

  par(mfrow = c(4, 3),
      mar = c(1.2, 1.5, 0.2, 0.2),
      oma = c(0.5, 3, 3.5, 0.5))

  for (i in 1:(max_delay)) {
    print(sprintf("dim(RTCC)==%s, N_days=%s, i=%s", dim(diagData$rep_tri_conf_cases), N_days, i))
    print(dim(diagData$rep_tri_conf_cases))
    plot(1:(N_days - i + 1),
         otcm_mu[1:(N_days - i + 1), i],
         ylim = c(0,
                  max(c(otcm_mu[1:(N_days - i + 1), i],
                        diagData$rep_tri_conf_cases[1:(N_days - i + 1), i]))) * 1.12,
         xlim = c(diagData$N_days_before, datList$N_days),
         type = "l",
         col = cls[i],
         axes = F,
         xlab = "",
         ylab = "")

    axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
    axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))
    box()

    text(diagData$N_days_before,
         max(c(otcm_mu[1:(N_days - i + 1), i],
               diagData$rep_tri_conf_cases[1:(N_days - i + 1), i])) * 1.06,
         paste0(i - 1, " days delay"),
         pos = 4,
         offset = -0.2,
         font = 2,
         cex = 1.15)

    lines(1:(N_days - i + 1),
          diagData$rep_tri_conf_cases[1:(N_days - i + 1), i],
          col = cls[i],
          lty = "12")
    points(1:(N_days - i + 1),
           diagData$rep_tri_conf_cases[1:(N_days - i + 1), i],
           col = cls[i],
           pch = 16,
           cex = 0.6)
  }
  mtext("Case count (N)", 2, 1, outer = T)
  mtext("Fit to data: Case time-series, stratified by days delay (points = data, lines = model)",
        3, 1, outer = T)
  dev.off()
}

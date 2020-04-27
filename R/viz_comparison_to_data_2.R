#' Comparison to data figure #2
#'
#' @inheritParams viz_comparison_to_data_1
#'
#' @return Side effects of plotting
#' @export
viz_comparison_to_data_2 <- function(cc) {
  samps    <- cc$extract
  datList  <- cc$config
  diagData <- cc$HMMMMMMMM # need to fix this...

  ########### Comparison to data 2 ###########

  par(mfrow = c(1, 1), mar = c(2.5, 3.5, 2, 0.2), oma = c(0, 0, 0, 0))
  N_days <- datList$N_days
  N_days_tot <- datList[["N_days"]] + datList[["N_days_extra"]]
  max_delay <- max(diagData$days_delay) + 1

  otcm_rep_del <- samps[["rep_tri_conf_cases_mu"]]
  for (k in 1:dim(otcm_rep_del)[1]) {
    for (i in 1:N_days) {
      for (j in 1:max_delay) {
        if ((i + j) >= (N_days + 2)) {
          otcm_rep_del[k, i, j] = 0
        }
      }
    }
  }

  otcm_interval <- apply(apply(samps[["rep_tri_conf_cases_mu"]], 1:2, sum), 2,
                         function(x) quantile(x, c(1, 39)/40))
  otcm_interval2 <- apply(apply(otcm_rep_del, 1:2, sum), 2, function(x)
                          quantile(x, c(1, 39)/40))
  otcm_mu <- apply(samps[["rep_tri_conf_cases_mu"]], 2:3, mean)
  otcm_mu2 <- apply(otcm_rep_del, 2:3, mean)

  plot(1:N_days, rowSums(otcm_mu), axes = F, xlab = "", ylab = "", type = "l",
       lty = "11", col = NA, ylim = c(0, max(otcm_interval)),
       xlim = c(which(rowSums(diagData$rep_tri_conf_cases) > 0)[1], N_days))

  polygon(c(1:N_days, N_days:1),
          c(otcm_interval[1, ],
            otcm_interval[2, N_days:1]),
          border = mTrsp("red", 30),
          col = mTrsp("red", 60))

  polygon(c(1:N_days, N_days:1),
          c(otcm_interval2[1, ],
            otcm_interval2[2, N_days:1]),
          border = mTrsp("navy", 30),
          col = mTrsp("navy", 60))

  lines(1:N_days, rowSums(otcm_mu), col = "red")
  lines(1:N_days, rowSums(otcm_mu2), col = "navy")
  lines(1:N_days, rowSums(diagData$rep_tri_conf_cases), col = mTrsp("navy", 150), lty = "12")

  points(1:N_days, rowSums(diagData$rep_tri_conf_cases), col = "navy", pch = 16, cex = 0.45)
  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))

  box()

  mtext("Case count (N)", 2, 2.3)
  mtext("Days since start", 1, 1.3)
  mtext("Fit to data: Cases diagnosed per day", 3, 0.6)
  legend("topleft",
         c("Model: new diagnoses (all)",
           "Model: new diagnoses (reported)",
           "Data: new diagnoses (reported)"), 
         col = c("red", "navy", "navy"),
         lty = c(1, 1, 3),
         pch = c(NA, NA, 16),
         pt.cex = 0.8,
         seg.len = 1.5,
         lwd = c(2, 1, 1),
         cex = 0.9,
         bty = "n",
         ncol = 1)

  new_inf <- apply(samps[["new_inf"]], 2, mean)
  new_sym <- apply(samps[["new_sym"]], 2, mean)
  new_hos <- apply(samps[["new_hos"]], 2, mean)
  new_die <- apply(samps[["new_die"]], 2, mean)
  new_res <- apply(samps[["res_inf"]] + samps[["res_sym"]] + samps[["res_hos"]], 2, mean)

  mx <- max(new_inf, new_sym, new_die, new_res)

  plot(1:N_days,
       rowSums(otcm_mu),
       axes = F, xlab = "", ylab = "",
       type = "l", lty = "11", col = NA,
       ylim = c(1, mx * 1.1),
       xlim = c(which(rowSums(diagData$rep_tri_conf_cases) > 0)[1], N_days_tot),
       log = "y")

  lines(1:N_days_tot, new_inf, col = "grey40", lwd = 2)
  lines(1:N_days_tot, new_sym, col = 2, lwd = 2)
  lines(1:N_days_tot, new_hos, col = "forestgreen", lwd = 2)
  lines(1:N_days_tot, new_die, col = "purple", lwd = 2)
  lines(1:N_days_tot, new_res, col = "brown", lwd = 2)
  lines(1:N_days, rowSums(otcm_mu), col = "orange", lwd = 2)
  lines(1:N_days, rowSums(otcm_mu2), col = "navy", lwd = 2)
  lines(1:N_days, rowSums(diagData$rep_tri_conf_cases), col = mTrsp("navy", 150), lty = "12")

  points(1:N_days, rowSums(diagData$rep_tri_conf_cases), col = "navy", pch = 16, cex = 0.6)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))
  box()

  mtext("Case count (N)", 2, 2.3)
  mtext("Days since start", 1, 1.3)
  mtext("Fit to data: Incident outcomes compared to data", 3, 0.6)

  legend("topleft",
         c("Model: new infections",
           "Model: new symptomatics",
           "Model: new hospitalizations",
           "Model: deaths",
           "Model: recoveries",
           "Model: new diagnoses (all)",
           "Model: new diagnoses (reported)",
           "Data: new diagnoses (reported)"), 
         col = c("grey40", "red", "forestgreen", "blue",
                 "purple", "brown", "orange", "navy", "navy"),
         lty = c(rep(1, 8), NA),
         pch = c(rep(NA, 8), 16),
         pt.cex = 0.7,
         seg.len = 1.2,
         lwd = 2,
         cex = 0.9,
         bty = "n",
         ncol = 2)
}

comparison_to_data_2 <- function(samps, datList) {
  ########### Comparison to data 2 ###########

  par(mfrow = c(1, 1), mar = c(2.5, 3.5, 2, 0.2), oma = c(0, 0, 0, 0))

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
       xlim = c(which(rowSums(rep_tri_conf_cases) > 0)[1], N_days))

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
  lines(1:N_days, rowSums(rep_tri_conf_cases), col = mTrsp("navy", 150), lty = "12")

  points(1:N_days, rowSums(rep_tri_conf_cases), col = "navy", pch = 16, cex = 0.45)
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
  ##### 

}

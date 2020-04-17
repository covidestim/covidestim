#' Current and cumulative outcomes figure
#'
#' @inheritParams viz_comparison_to_data_1
#'
#' @return Side effects of plotting
#' @export
viz_current_and_cumulative_outcomes <- function(samps, datList, ...) {

  smp  <- ceiling(seq(1,nrow( samps[["new_inf"]]),length.out=10))
  N_days_tot <- datList[["N_days"]] + datList[["N_days_extra"]]

  ########### Current outcomes and cumulative final outcomes ###########
  par(mfrow = c(2, 3), mar = c(1.2, 2.1, 0.2, 0.2), oma = c(3.5, 3, 3.5, 0.5))

  #### Current Asymptomatic
  otcm <- samps[["cur_inf"]]
  plot(1:N_days_tot,
       apply(otcm, 2, mean),
       type = "l",
       ylim = c(0, quantile(otcm[, N_days_tot], 0.975)) * 1.04,
       las = 1, main = "", xlab = "", ylab = "", axes = F)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))
  box()

  polygon(c(1:N_days_tot, N_days_tot:1),
          c(apply(otcm, 2, function(x) quantile(x, 0.025)),
            apply(otcm, 2, function(x) quantile(x, 0.975))[N_days_tot:1]),
          border = "grey82", col = "grey90")

  for (i in smp)
    lines(1:N_days_tot, otcm[i, ], col = mTrsp(3, 120))

  lines(1:N_days_tot, apply(otcm, 2, mean))
  text(0,
       quantile(otcm[, N_days_tot], 0.975) * 1.01,
       "Current asymptomatic",
       pos = 4, offset = 0.2, font = 2, cex = 1.15)

  #### Current Symptomatic, Non-Hospitalized
  otcm <- samps[["cur_sym"]]

  plot(1:N_days_tot,
       apply(otcm, 2, mean),
       type = "l",
       ylim = c(0, quantile(otcm[, N_days_tot], 0.975)) * 1.04,
       las = 1, main = "", xlab = "", ylab = "", axes = F)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))

  box()

  polygon(c(1:N_days_tot, N_days_tot:1),
          c(apply(otcm, 2, function(x) quantile(x, 0.025)),
            apply(otcm, 2, function(x) quantile(x, 0.975))[N_days_tot:1]),
          border = "grey82", col = "grey90")

  for (i in smp)
    lines(1:N_days_tot, otcm[i, ], col = mTrsp(3, 120))

  lines(1:N_days_tot, apply(otcm, 2, mean))

  text(0,
       quantile(otcm[, N_days_tot], 0.975) * 1.01,
       "Current symptomatic, non-hospitalized", 
       pos = 4, offset = 0.2, font = 2, cex = 1.15)

  #### Current Hospitalized
  otcm <- samps[["cur_hos"]]

  plot(1:N_days_tot,
       apply(otcm, 2, mean),
       type = "l",
       ylim = c(0, quantile(otcm[, N_days_tot], 0.975)) * 1.04,
       las = 1, main = "", xlab = "", ylab = "", axes = F)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))

  box()

  polygon(c(1:N_days_tot, N_days_tot:1),
          c(apply(otcm, 2, function(x) quantile(x, 0.025)),
            apply(otcm, 2, function(x) quantile(x, 0.975))[N_days_tot:1]),
          border = "grey82", col = "grey90")

  for (i in smp)
    lines(1:N_days_tot, otcm[i, ], col = mTrsp(3, 120))

  lines(1:N_days_tot, apply(otcm, 2, mean))

  text(0,
       quantile(otcm[, N_days_tot], 0.975) * 1.01,
       "Current hospitalized",
       pos = 4, offset = 0.2, font = 2, cex = 1.15)

  #### Cumulative dead
  otcm <- samps[["cum_die"]]
  plot(1:N_days_tot,
       apply(otcm, 2, mean),
       type = "l",
       ylim = c(0, quantile(otcm[, N_days_tot], 0.975)) * 1.04,
       las = 1, main = "", xlab = "", ylab = "", axes = F)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))
  box()

  polygon(c(1:N_days_tot, N_days_tot:1),
          c(apply(otcm, 2, function(x) quantile(x, 0.025)),
            apply(otcm, 2, function(x) quantile(x, 0.975))[N_days_tot:1]),
          border = "grey82", col = "grey90")

  for (i in smp)
    lines(1:N_days_tot, otcm[i, ], col = mTrsp(3, 120))

  lines(1:N_days_tot, apply(otcm, 2, mean))

  text(0,
       quantile(otcm[, N_days_tot], 0.975) * 1.01,
       "Cumulative dead", pos = 4, offset = 0.2, font = 2, cex = 1.15)

  #### Cumulative recovered
  otcm <- samps[["cum_res"]]

  plot(1:N_days_tot,
       apply(otcm, 2, mean),
       type = "l",
       ylim = c(0, quantile(otcm[, N_days_tot], 0.975)) * 1.04,
       las = 1, main = "", xlab = "", ylab = "", axes = F)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))
  box()

  polygon(c(1:N_days_tot, N_days_tot:1),
          c(apply(otcm, 2, function(x) quantile(x, 0.025)),
            apply(otcm, 2, function(x) quantile(x, 0.975))[N_days_tot:1]),
          border = "grey82", col = "grey90")

  for (i in smp)
    lines(1:N_days_tot, otcm[i, ], col = mTrsp(3, 120))

  lines(1:N_days_tot, apply(otcm, 2, mean))

  text(0,
       quantile(otcm[, N_days_tot], 0.975) * 1.01,
       "Cumulative recovered",
       pos = 4, offset = 0.2, font = 2, cex = 1.15)

  mtext("Count of outcomes (N)", 2, 1, outer = T)
  mtext("Current and cumulative outcomes: new infections, symptomatics, hospitalizations, deaths, and recoveries",
        3, 1, outer = T)
  mtext("Days since start", 1, 1, outer = T)
}

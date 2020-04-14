priors_vs_posteriors <- function(samps, datList) {
  ############ Priors vs. Posteriors ###########

  par(mfrow = c(3, 4),
      mar = c(1.2, 1.3, 2, 0.2),
      oma = c(0.5, 2.5, 0.5, 0.5))

  ## log_new_inf_0
  mu = datList[["pri_log_new_inf_0_mu"]]
  sd = datList[["pri_log_new_inf_0_sd"]]

  hist(samps[["log_new_inf_0"]],
       col = "goldenrod",
       border = F,
       breaks = 20,
       axes = F,
       main = "",
       xlab = "",
       ylab = "",
       xlim = range(c(samps[["log_new_inf_0"]], mu - sd/2, mu + sd/2)), probability = T)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))
  box()

  zz <- range(c(samps[["log_new_inf_0"]], mu - sd, mu + sd))
  x_pts <- seq(zz[1], zz[2], length.out = 101)
  y_pts <- dnorm(x_pts, mu, sd)
  lines(x_pts, y_pts, col = mTrsp(4, 200), lwd = 2)
  mtext("log_new_inf_0", 3, 0.5, font = 2, cex = 0.8)

  ## log_new_inf_drift

  mu = datList[["pri_log_new_inf_drift_mu"]]
  sd = datList[["pri_log_new_inf_drift_sd"]]
  mod = samps[["log_new_inf_drift"]]
  hist(mod,
       col = "goldenrod",
       border = F,
       breaks = 10,
       axes = F,
       main = "",
       xlab = "",
       ylab = "",
       xlim = range(c(mod, mu - sd/2, mu + sd/2)),
       probability = T)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))

  box()

  zz <- range(c(mod, mu - sd, mu + sd))
  x_pts <- seq(zz[1], zz[2], length.out = 101)
  y_pts <- dnorm(x_pts, mu, sd)

  lines(x_pts, y_pts, col = mTrsp(4, 200), lwd = 2)
  mtext("log_new_inf_drift", 3, 0.5, font = 2, cex = 0.8)
  ## pri_sigma_deriv1_log_new_inf_sd
  mu = 0
  sd = datList[["pri_sigma_deriv1_log_new_inf_sd"]]
  mod = samps[["sigma_deriv1_log_new_inf"]]
  hist(mod,
       col = "goldenrod",
       border = F,
       breaks = 20,
       axes = F, main = "", xlab = "", ylab = "",
       xlim = range(c(mod, 0, mu + sd/2)),
       probability = T)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))
  box()

  zz <- range(c(mod, mu - sd, mu + sd))
  x_pts <- seq(zz[1], zz[2], length.out = 101)
  y_pts <- dnorm(x_pts, mu, sd)

  lines(x_pts, y_pts, col = mTrsp(4, 200), lwd = 2)
  mtext("sigma_deriv1_log_new_inf", 3, 0.5, font = 2, cex = 0.8)

  ## p_sym_if_inf
  mn = datList[["pri_p_sym_if_inf_mn"]]
  ss = datList[["pri_p_sym_if_inf_ss"]]
  mod = samps[["p_sym_if_inf"]]
  hist(mod,
       col = "goldenrod",
       border = F,
       breaks = 15,
       axes = F,
       main = "",
       xlab = "",
       ylab = "",
       xlim = range(c(mod, qbeta(c(2e-04, 0.9998), mn * ss, (1 - mn) * ss))),
       probability = T)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))

  box()
  zz <- range(c(mod, 0, 1))
  x_pts <- seq(zz[1], zz[2], length.out = 101)
  y_pts <- dbeta(x_pts, mn * ss, (1 - mn) * ss)
  lines(x_pts, y_pts, col = mTrsp(4, 200), lwd = 2)
  mtext("p_sym_if_inf", 3, 0.5, font = 2, cex = 0.8)

  ## p_hos_if_sym
  mn = datList[["pri_p_hos_if_sym_mn"]]
  ss = datList[["pri_p_hos_if_sym_ss"]]
  mod = samps[["p_hos_if_sym"]]
  hist(mod,
       col = "goldenrod",
       border = F,
       breaks = 15,
       axes = F,
       main = "",
       xlab = "",
       ylab = "",
       xlim = range(c(mod, qbeta(c(2e-04, 0.9998), mn * ss, (1 - mn) * ss))),
       probability = T)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))

  box()
  zz <- range(c(mod, 0, 1))
  x_pts <- seq(zz[1], zz[2], length.out = 101)
  y_pts <- dbeta(x_pts, mn * ss, (1 - mn) * ss)
  lines(x_pts, y_pts, col = mTrsp(4, 200), lwd = 2)
  mtext("p_hos_if_sym", 3, 0.5, font = 2, cex = 0.8)

  ## p_die_if_hos
  mn = datList[["pri_p_die_if_hos_mn"]]
  ss = datList[["pri_p_die_if_hos_ss"]]
  mod = samps[["p_die_if_hos"]]
  hist(mod,
       col = "goldenrod",
       border = F,
       breaks = 15,
       axes = F,
       main = "",
       xlab = "",
       ylab = "",
       xlim = range(c(mod, qbeta(c(2e-04, 0.9998), mn * ss, (1 - mn) * ss))),
       probability = T)

  axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
  axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))

  box()
  zz <- range(c(mod, 0, 1))
  x_pts <- seq(zz[1], zz[2], length.out = 101)
  y_pts <- dbeta(x_pts, mn * ss, (1 - mn) * ss)
  lines(x_pts, y_pts, col = mTrsp(4, 200), lwd = 2)
  mtext("p_die_if_hos", 3, 0.5, font = 2, cex = 0.8)

  mtext("Model parameters: prior (blue) versus posterior (yellow)", 3, 1, outer = T)
}

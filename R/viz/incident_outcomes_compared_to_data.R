smp <- ceiling(seq(1, nrow(samps[["new_inf"]]), length.out = 10))

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
     xlim = c(which(rowSums(rep_tri_conf_cases) > 0)[1], N_days_tot),
     log = "y")

lines(1:N_days_tot, new_inf, col = "grey40", lwd = 2)
lines(1:N_days_tot, new_sym, col = 2, lwd = 2)
lines(1:N_days_tot, new_hos, col = "forestgreen", lwd = 2)
lines(1:N_days_tot, new_die, col = "purple", lwd = 2)
lines(1:N_days_tot, new_res, col = "brown", lwd = 2)
lines(1:N_days, rowSums(otcm_mu), col = "orange", lwd = 2)
lines(1:N_days, rowSums(otcm_mu2), col = "navy", lwd = 2)
lines(1:N_days, rowSums(rep_tri_conf_cases), col = mTrsp("navy", 150), lty = "12")

points(1:N_days, rowSums(rep_tri_conf_cases), col = "navy", pch = 16, cex = 0.6)

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

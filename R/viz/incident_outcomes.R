########### Incident Outcomes ########### New Infections
par(mfrow = c(2, 3),
    mar = c(1.2, 2, 0.2, 0.2),
    oma = c(3.5, 3, 3.5, 0.5))

otcm <- samps[["new_inf"]]

plot(1:N_days_tot,
     apply(otcm, 2, mean),
     type = "l",
     ylim = c(0, quantile(otcm[, N_days_tot], 0.975)) * 1.04,
     las = 1,
     main = "",
     axes = F,
     xlab = "",
     ylab = "")

axis(2, las = 1, tcl = -0.07, mgp = c(3, 0.15, 0))
axis(1, las = 1, tcl = -0.07, mgp = c(3, 0.1, 0))

box()

polygon(c(1:N_days_tot, N_days_tot:1),
        c(apply(otcm, 2, function(x) quantile(x, 0.025)),
          apply(otcm, 2, function(x) quantile(x, 0.975))[N_days_tot:1]),
        border = "grey82",
        col = "grey90")

for (i in smp)
  lines(1:N_days_tot, otcm[i, ], col = mTrsp(3, 120))

lines(1:N_days_tot, apply(otcm, 2, mean))

text(0,
     quantile(otcm[, N_days_tot], 0.975) * 1.01,
     "New infections per day",
     pos = 4, offset = 0.2, font = 2, cex = 1.15)

#### New Symp
otcm <- samps[["new_sym"]]

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

text(0, quantile(otcm[, N_days_tot], 0.975) * 1.01,
     "New symptomatic per day",
     pos = 4, offset = 0.2, font = 2, cex = 1.15)

#### New Hosp
otcm <- samps[["new_hos"]]
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
        border = "grey82",
        col = "grey90")

for (i in smp)
  lines(1:N_days_tot, otcm[i, ], col = mTrsp(3, 120))

lines(1:N_days_tot, apply(otcm, 2, mean))

text(0,
     quantile(otcm[, N_days_tot], 0.975) * 1.01,
     "New hospitalizations per day",
     pos = 4, offset = 0.2, font = 2, cex = 1.15)

#### New death
otcm <- samps[["new_die"]]
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
     "Deaths per day",
     pos = 4, offset = 0.2, font = 2, cex = 1.15)

#### New recovery
otcm <- samps[["res_inf"]] + samps[["res_sym"]] + samps[["res_hos"]]

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

text(0, quantile(otcm[, N_days_tot], 0.975) * 1.01,
     "Recoveries per day",
     pos = 4, offset = 0.2, font = 2, cex = 1.15)

mtext("Count of outcomes (N)", 2, 1, outer = T)
mtext("Incident outcomes: new infections, symptomatics, hospitalizations, deaths, and recoveries",
      3, 1, outer = T)
mtext("Days since start", 1, 1, outer = T)



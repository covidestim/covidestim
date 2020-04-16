  # Delay from infection to symptoms
  # https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported
  # Median incubation 5.1 days (4.5, 5.8 );
  # 97.5 percentile 11.5 days (8.2, 15.6);
  # 2.5 percentile 2.2 days (1.8, 2.9)

  # functions to get parameters from the Lauer paper
  exp(
    stats::optim(
      c(-1, -1),
      function(zz) {
        z   <- exp(zz)
        prd <- qgamma(c(1, 20, 39)/40, 1/z[1]/z[1], 1/z[1]/z[1]/z[2])
        sum( (prd - c(2.2, 5.1, 11.5))^2 * c(1, 5, 1) )
      },
      method="BFGS"
    )$par
  ) -> opt_par1

  # qgamma(c(1, 20, 39)/40, opt_par1[1]^-2, opt_par1[1]^-2/opt_par1[2])

  stats::optimize(
    function(z) {
      prd <- qgamma(c(1, 39)/40, 1/z[1]/z[1], 1/z[1]/z[1]/5.1)
      sum((prd - c(4.5, 5.8))^2)
    },
    c(0.01, 1)
  )$minimum -> opt_par2

  # qgamma(c(1, 39)/40, opt_par2^-2, opt_par2^-2/5.1)

  gammapar <- function(tgt) {
     tgt <- as.numeric(tgt)
     mn <- tgt[1]
     cir <- (tgt[3] - tgt[2])
     xopt <- function(b, mn, cir) {
       cir2 <- qgamma(c(1, 39)/40, mn * b, b)
       cir2 <- cir2[2] - cir2[1]
       (cir2 - cir)^2
     }
     zz <- optimize(xopt, c(0.1, 1e+05), mn = mn, cir = cir)$minimum
     c(zz * mn, zz)
   }
   
   
   

# tests

pars_of_interest <- c("p_die_if_sym", "p_die_if_sev", "p_sym_if_inf", 
                      "p_sev_if_sym","p_diag_if_sym", "p_diag_if_sev", 
                      "scale_dx_delay_sym", "scale_dx_delay_sev", 
                      "log_new_inf_0")

rstan::check_hmc_diagnostics(ccr$result)
summary(as_tibble(ccr$summary)$Rhat)
summary(as_tibble(ccr$summary)$n_eff)

rstan::stan_diag(ccr$result,
                 information = c("sample","stepsize", "treedepth","divergence"),
                 chain = 0)
rstan::traceplot(ccr$result, pars = pars_of_interest)

samp <- ccr$extracted
pairs(cbind(samp$p_die_if_sym, samp$p_die_if_sev, 
            samp$p_sym_if_inf, samp$p_sev_if_sym, 
            samp$p_diag_if_sym, samp$p_diag_if_sev, 
            samp$scale_dx_delay_sym, samp$scale_dx_delay_sev, 
            samp$log_new_inf_0),
      labels = pars_of_interest, 
      pch=16, col=4, cex=.3)
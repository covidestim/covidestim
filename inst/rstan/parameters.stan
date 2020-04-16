parameters {
  
// INCIDENCE  

//  real                  log_new_inf_0; // intercept in log space
//  real                  log_new_inf_drift; // mean day on day change
//  vector[N_days_tot-1]  deriv1_log_new_inf; // first derivative of the random walk 
//  real<lower=0>         sigma_deriv1_log_new_inf; // parameter for the SD of the rw

  vector[n_spl_par]     b_spline;

// SYMPTOMS AND CARE 
  real<lower=0, upper=1>  p_sym_if_inf;
  real<lower=0, upper=1>  p_hos_if_sym;
  real<lower=0, upper=1>  p_die_if_hos;
// delay associated with transition to next illness state
  real<lower=0>           inf_prg_delay_mn;
  real<lower=0>           sym_prg_delay_mn;
  real<lower=0>           hos_prg_delay_mn;
// delay associated with transition to recovered 
  real<lower=0>           inf_res_delay_mn;
  real<lower=0>           sym_res_delay_mn;
  real<lower=0>           hos_res_delay_mn;
  
// DIAGNOSIS // probability of diagnosis at each illness state
  real<lower=0, upper=1>  p_diag_if_inf;
  real<lower=0, upper=1>  p_diag_if_sym;
  real<lower=0, upper=1>  p_diag_if_hos;
  
// REPORTING DElAYS  
  real<lower=0>           report_delay_mn;
  
// LIKELIHOOD
//  real<lower=0>           inv_sqrt_phi;
}

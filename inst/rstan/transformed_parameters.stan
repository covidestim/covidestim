transformed parameters {
///~~~~~~~ Define ~~~~~~~
// INCIDENCE
  vector[N_days_tot]      log_new_inf;
  vector[N_days_tot]      new_inf;

  vector[N_days_tot-2]   deriv2_log_new_inf;

  //vector[n_spl_par-2]     deriv2_b_spline;
  //vector[n_spl_par-1]     deriv1_b_spline;
 
// DELAYS probabilitily of reporting w x days delay (PDF of delay distribution)
  vector[N_days_tot]  inf_prg_delay;
  vector[N_days_tot]  sym_prg_delay;
  vector[N_days_tot]  hos_prg_delay;
  vector[N_days_tot]  inf_res_delay;
  vector[N_days_tot]  sym_res_delay;
  vector[N_days_tot]  hos_res_delay;
  vector[N_days_tot]  report_delay;

// OUTCOMES
  vector[N_days_tot]  new_sym; // new in state per day
  vector[N_days_tot]  new_hos;
  vector[N_days_tot]  new_die;
  vector[N_days_tot]  res_inf; // resolve from state (all back to no disease)
  vector[N_days_tot]  res_sym;
  vector[N_days_tot]  res_hos;
  vector[N_days_tot]  res_all;
  
  vector[N_days_tot]  new_sym_u; // enter new state << undiagnosed >>
  vector[N_days_tot]  new_hos_u;
  vector[N_days_tot]  new_die_u;
  vector[N_days_tot]  res_inf_u; // resolve from state << undiagnosed >>
  vector[N_days_tot]  res_sym_u;
  vector[N_days_tot]  res_hos_u;
  
  vector[N_days_tot]  cum_inf; // cumulative in state
  vector[N_days_tot]  cum_sym; 
  vector[N_days_tot]  cum_hos;
  vector[N_days_tot]  cum_die;
  vector[N_days_tot]  cum_res;
  vector[N_days_tot]  cur_inf; // currently in state 
  vector[N_days_tot]  cur_sym;
  vector[N_days_tot]  cur_hos;
  
  vector[N_days_tot]  cum_diag_inf; // cumulative diagnoses by state
  vector[N_days_tot]  cum_diag_sym;
  vector[N_days_tot]  cum_diag_hos;
  vector[N_days_tot]  cur_inf_u; // currently in state << undiagnosed >>
  vector[N_days_tot]  cur_sym_u;
  vector[N_days_tot]  cur_hos_u;

  vector[N_days_tot]  diag_inf; // incident diagnosis
  vector[N_days_tot]  diag_sym;
  vector[N_days_tot]  diag_hos;
  vector[N_days_tot]  diag_all;
  
//  real                phi;
  matrix[N_days,Max_delay+1]  rep_tri_conf_cases_mu;
  
///~~~~~~~ Assign values ~~~~~~~
//  phi = pow(inv_sqrt_phi,-2);
  
// INCIDENCE

//if(rw_yes == 1){
   log_new_inf[1] = log_new_inf_0;
   for(i in 1:(N_days_tot-1)) {
    log_new_inf[i+1] =  log_new_inf[i] + deriv1_log_new_inf[i] ;
  }
  new_inf = exp(log_new_inf);
  for(i in 1:(N_days_tot-2)) {
    deriv2_log_new_inf[i] = log_new_inf[i+1] * 2 - log_new_inf[i] - log_new_inf[i+2];
  }
//} else {
 // log_new_inf = spl_basis * b_spline;
//  new_inf  = exp(log_new_inf);
// second derivative
//for(i in 1:(n_spl_par-2)) {
//  deriv2_b_spline[i] = b_spline[i+1] * 2 - b_spline[i] - b_spline[i+2];
//}
// first derivative
//for(i in 1:(n_spl_par-1)) {
//  deriv1_b_spline[i] = b_spline[i+1] - b_spline[i];
//}
//}

// DELAYS // progression
  for(i in 1:N_days_tot) {
    inf_prg_delay[i] = gamma_cdf(i+0.0, inf_prg_delay_shap, inf_prg_delay_rate)
      - gamma_cdf(i-1.0, inf_prg_delay_shap, inf_prg_delay_rate);
    sym_prg_delay[i] = gamma_cdf(i+0.0, sym_prg_delay_shap, sym_prg_delay_rate)
      - gamma_cdf(i-1.0, sym_prg_delay_shap, sym_prg_delay_rate);
    hos_prg_delay[i] = gamma_cdf(i+0.0, hos_prg_delay_shap, hos_prg_delay_rate)
      - gamma_cdf(i-1.0, hos_prg_delay_shap, hos_prg_delay_rate);
  }
  
// DELAYS // resolution of case
  for(i in 1:N_days_tot) {
    inf_res_delay[i] = gamma_cdf(i+0.0, inf_res_delay_shap, inf_res_delay_rate)
      - gamma_cdf(i-1.0, inf_res_delay_shap, inf_res_delay_rate);
    sym_res_delay[i] = gamma_cdf(i+0.0, sym_res_delay_shap, sym_res_delay_rate)
      - gamma_cdf(i-1.0, sym_res_delay_shap, sym_res_delay_rate);
    hos_res_delay[i] = gamma_cdf(i+0.0, hos_res_delay_shap, hos_res_delay_rate)
      - gamma_cdf(i-1.0, hos_res_delay_shap, hos_res_delay_rate);
  }
  
// DELAYS // reporting // needs to be simplified 
  for(i in 1:N_days_tot) {
    report_delay[i] = gamma_cdf(i+0.0, pri_report_delay_shap, pri_report_delay_rate) -
      gamma_cdf(i-1.0, pri_report_delay_shap, pri_report_delay_rate);
  }
  report_delay = report_delay/sum(report_delay);
  
// CASCADE OF INCIDENT OUTCOMES (TOTAL)
  new_sym = rep_vector(0, N_days_tot);
  new_hos = rep_vector(0, N_days_tot);
  new_die = rep_vector(0, N_days_tot);
  res_inf = rep_vector(0, N_days_tot);
  res_sym = rep_vector(0, N_days_tot);
  res_hos = rep_vector(0, N_days_tot);
  
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        res_inf[i+(j-1)] += new_inf[i] * 
        (1-p_sym_if_inf) * inf_res_delay[j]; // probability resolved on day J
      }
    }
  }
  // could be 'prof_inf' <- infections that progress
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_sym[i+(j-1)] += new_inf[i] * p_sym_if_inf * inf_prg_delay[j];
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        res_sym[i+(j-1)] += new_sym[i] * (1-p_hos_if_sym) * sym_res_delay[j];
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_hos[i+(j-1)] += new_sym[i] * p_hos_if_sym * sym_prg_delay[j];
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        res_hos[i+(j-1)] += new_hos[i] * (1-p_die_if_hos) * hos_res_delay[j];
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_die[i+(j-1)] += new_hos[i] * p_die_if_hos * hos_prg_delay[j];
      }
    }
  }
  res_all = res_inf + res_sym + res_hos;
  
// CUMULATIVE
  cum_inf = cumulative_sum(new_inf);
  cum_sym = cumulative_sum(new_sym);
  cum_hos = cumulative_sum(new_hos);
  cum_die = cumulative_sum(new_die);
  cum_res = cumulative_sum(res_all);
  
// CURRENT STATUS // n in each state at a given time i 
  cur_inf = rep_vector(0, N_days_tot);
  cur_sym = rep_vector(0, N_days_tot);
  cur_hos = rep_vector(0, N_days_tot);
  cur_inf[1] = new_inf[1] - new_sym[1] - res_inf[1];
  for(i in 2:N_days_tot) {
    cur_inf[i] = cur_inf[i-1] + new_inf[i] - new_sym[i] - res_inf[i];
  }
  cur_sym[1] = new_sym[1] - new_hos[1] - res_sym[1];
  for(i in 2:N_days_tot) {
    cur_sym[i] = cur_sym[i-1] + new_sym[i] - new_hos[i] - res_sym[i];
  }
  cur_hos[1] = new_hos[1] - new_die[1] - res_hos[1];
  for(i in 2:N_days_tot) {
    cur_hos[i] = cur_hos[i-1] + new_hos[i] - new_die[i] - res_hos[i];
  }
  
// CASCADE OF INCIDENT OUTCOMES << UNDIAGNOSED >>
// checks after: make sure everything adds
  new_sym_u = rep_vector(0, N_days_tot);
  new_hos_u = rep_vector(0, N_days_tot);
  new_die_u = rep_vector(0, N_days_tot);
  res_inf_u = rep_vector(0, N_days_tot);
  res_sym_u = rep_vector(0, N_days_tot);
  res_hos_u = rep_vector(0, N_days_tot);
  
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        res_inf_u[i+(j-1)] += new_inf[i] * (1-p_sym_if_inf) * inf_res_delay[j]
          * pow(1-p_diag_if_inf,j-1);
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_sym_u[i+(j-1)] += new_inf[i] * p_sym_if_inf * inf_prg_delay[j] *
          pow(1-p_diag_if_inf,j-1);
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        res_sym_u[i+(j-1)] += new_sym_u[i] * (1-p_hos_if_sym) *
          sym_res_delay[j] * pow(1-p_diag_if_sym,j-1);
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_hos_u[i+(j-1)] += new_sym_u[i] * p_hos_if_sym * sym_prg_delay[j] *
          pow(1-p_diag_if_sym,j-1);
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        res_hos_u[i+(j-1)] += new_hos_u[i] * (1-p_die_if_hos) *
          hos_res_delay[j] * pow(1-p_diag_if_hos,j-1);
      }
    }
  }
  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_die_u[i+(j-1)] += new_hos_u[i] * p_die_if_hos * hos_prg_delay[j] *
          pow(1-p_diag_if_hos,j-1);
      }
    }
  }
  
// CURRENT STATUS  << UNDIAGNOSED >>
  cur_inf_u = rep_vector(0, N_days_tot);
  cur_sym_u = rep_vector(0, N_days_tot);
  cur_hos_u = rep_vector(0, N_days_tot);
  
  cur_inf_u[1] = new_inf[1] - new_sym_u[1] - res_inf_u[1];
  for(i in 2:N_days_tot) {
    cur_inf_u[i] = cur_inf_u[i-1]*(1-p_diag_if_inf) + new_inf[i] - new_sym_u[i]
      - res_inf_u[i];
  }
  cur_sym_u[1] = new_sym_u[1] - new_hos_u[1] - res_sym_u[1];
  for(i in 2:N_days_tot) {
    cur_sym_u[i] = cur_sym_u[i-1]*(1-p_diag_if_sym) + new_sym_u[i] -
      new_hos_u[i] - res_sym_u[i]; 
  }
  cur_hos_u[1] = new_hos_u[1] - new_die_u[1] - res_hos_u[1];  
  for(i in 2:N_days_tot) {
    cur_hos_u[i] = cur_hos_u[i-1]*(1-p_diag_if_hos) + new_hos_u[i] -
      new_die_u[i] - res_hos_u[i];
  }
  
// DIAGNOSIS // check diag_inf[i] = cur_inf[i] - cur_inf_u[i]
  for(i in 1:N_days_tot) {
    diag_inf[i] = cur_inf_u[i]*p_diag_if_inf;
    diag_sym[i] = cur_sym_u[i]*p_diag_if_sym;
    diag_hos[i] = cur_hos_u[i]*p_diag_if_hos;
  }
  diag_all = diag_inf + diag_sym + diag_hos;
// REPORTING TIRANGLE
  for(i in 1:N_days) {
    for(j in 1:(Max_delay+1)) {
      rep_tri_conf_cases_mu[i,j] = diag_all[i] * report_delay[j];
    }
  }
// CUMULATIVE DIAGNOSES from each state
  cum_diag_inf = cumulative_sum(diag_inf);
  cum_diag_sym = cumulative_sum(diag_sym);
  cum_diag_hos = cumulative_sum(diag_hos);
}


data {
  ///~~~~~~~ Define ~~~~~~~
  // INPUT DATA
  int<lower=0>           N_days;
  int<lower=0>           N_days_before;
  int<lower=0>           Max_delay;
  int<lower=0>           N_spl_par;
  int<lower=0>           obs_cas[N_days]; 
  int<lower=0>           obs_die[N_days];
  matrix[N_days+N_days_before,N_spl_par] spl_basis;
  
  int<lower=0>           N_spl_par_dx;
  matrix[N_days+N_days_before,N_spl_par_dx] spl_basis_dx;

  // PRIORS
  // incidence time-series 
  real                   pri_log_new_inf_0_mu;
  real<lower=0>          pri_log_new_inf_0_sd;
  real                   pri_logRt_mu;
  real<lower=0>          pri_logRt_sd;
  real                   pri_serial_i_a;
  real<lower=0>          pri_serial_i_b;
  real                   pri_inf_imported_mu;
  real<lower=0>          pri_inf_imported_sd;
  real<lower=0>          pri_deriv1_spl_par_sd;
  real<lower=0>          pri_deriv2_spl_par_sd;
  
  // p(progression)
  real<lower=0>          pri_p_sym_if_inf_a; 
  real<lower=0>          pri_p_sym_if_inf_b;
  real<lower=0>          pri_p_sev_if_sym_a;
  real<lower=0>          pri_p_sev_if_sym_b;
  real<lower=0>          pri_p_die_if_sev_a;
  real<lower=0>          pri_p_die_if_sev_b;
  real<lower=0>          pri_p_die_if_inf_a;
  real<lower=0>          pri_p_die_if_inf_b;
  // p(diag)  
  real<lower=0>          pri_rr_diag_sym_vs_sev_a;
  real<lower=0>          pri_rr_diag_sym_vs_sev_b;
  real<lower=0>          pri_p_diag_if_sev_a;
  real<lower=0>          pri_p_diag_if_sev_b;
  // delay to progression
  real<lower=0>          inf_prg_delay_shap;
  real<lower=0>          inf_prg_delay_rate;
  real<lower=0>          sym_prg_delay_shap;
  real<lower=0>          sym_prg_delay_rate;
  real<lower=0>          sev_prg_delay_shap;
  real<lower=0>          sev_prg_delay_rate;
  //
  real<lower=0>          scale_dx_delay_sym_a; 
  real<lower=0>          scale_dx_delay_sym_b; 
  real<lower=0>          scale_dx_delay_sev_a; 
  real<lower=0>          scale_dx_delay_sev_b;
  
  // delay to report
  real<lower=0>          pri_cas_rep_delay_shap_a;
  real<lower=0>          pri_cas_rep_delay_shap_b;
  real<lower=0>          pri_cas_rep_delay_rate_a;
  real<lower=0>          pri_cas_rep_delay_rate_b;
  
  real<lower=0>          pri_die_rep_delay_shap_a;
  real<lower=0>          pri_die_rep_delay_shap_b;
  real<lower=0>          pri_die_rep_delay_rate_a;
  real<lower=0>          pri_die_rep_delay_rate_b;
  
  real<lower=0>          pri_inv_sqrt_phi;

  int<lower = 0, upper = 1> cas_yes; 
  int<lower = 0, upper = 1> die_yes; 
  int<lower = 0, upper = 1> obs_cas_rep; 
  int<lower = 0, upper = 1> obs_die_rep; 
  int<lower = 1, upper = 7> n_day_av; 
}
///////////////////////////////////////////////////////////
transformed data {
  int  N_days_tot;
  int  nda0;
  vector[Max_delay]  inf_prg_delay;
  vector[Max_delay]  sym_prg_delay;
  vector[Max_delay]  sev_prg_delay;
  vector[Max_delay]  inf_cum_prg_delay;
  vector[Max_delay]  sym_cum_prg_delay; 
  vector[Max_delay]  sev_cum_prg_delay; 
  
  N_days_tot = N_days + N_days_before; 
  nda0 = n_day_av - 1;
  
  //  progression
  for(i in 1:Max_delay) {
    inf_prg_delay[i] = gamma_cdf(i+0.0, inf_prg_delay_shap, inf_prg_delay_rate)
      - gamma_cdf(i-1.0, inf_prg_delay_shap, inf_prg_delay_rate);
    sym_prg_delay[i] = gamma_cdf(i+0.0, sym_prg_delay_shap, sym_prg_delay_rate)
      - gamma_cdf(i-1.0, sym_prg_delay_shap, sym_prg_delay_rate);
    sev_prg_delay[i] = gamma_cdf(i+0.0, sev_prg_delay_shap, sev_prg_delay_rate)
      - gamma_cdf(i-1.0, sev_prg_delay_shap, sev_prg_delay_rate);
  }
  
  inf_cum_prg_delay = cumulative_sum(inf_prg_delay);
  sym_cum_prg_delay = cumulative_sum(sym_prg_delay);
  sev_cum_prg_delay = cumulative_sum(sev_prg_delay);

}
///////////////////////////////////////////////////////////
parameters {
// INCIDENCE  
  real                    log_new_inf_0; // starting intercept
  real<lower=0>           serial_i; // serial interval
  real<lower=0>           inf_imported; // imported cases
  vector[N_spl_par]       spl_par;

// SYMPTOMS AND CARE 
  real<lower=0, upper=1>    p_sym_if_inf;
  real<lower=0, upper=1>  p_sev_if_sym;
  real<lower=0, upper=1>  p_die_if_sev;
  real<lower=0, upper=1>   scale_dx_delay_sym; 
  real<lower=0, upper=1>   scale_dx_delay_sev; 
  
  real<lower=0.5>     cas_rep_delay_shap;
  real<lower=0.05>    cas_rep_delay_rate;
  real<lower=0.5>     die_rep_delay_shap;
  real<lower=0.05>    die_rep_delay_rate;

// DIAGNOSIS // probability of diagnosis at each illness state
  real<lower=0, upper=1>  p_diag_if_sev;
  vector[N_spl_par_dx]    spl_par_dx;

// LIKELIHOOD
  real<lower=0>           inv_sqrt_phi_c;
  real<lower=0>           inv_sqrt_phi_d;

}
///////////////////////////////////////////
transformed parameters {
///~~~~~~~ Define ~~~~~~~
// INCIDENCE
  vector[N_days_tot]      log_new_inf;
  vector[N_days_tot]      new_inf;
  vector[N_days_tot]      deriv1_log_new_inf;
  //~~~~~~
  vector[N_days_tot]      logRt;
  vector[N_days_tot]      Rt;
  vector[N_spl_par-1]     deriv1_spl_par;
  vector[N_spl_par-2]     deriv2_spl_par;

  real  p_die_if_inf;
  vector[N_days_tot]  p_diag_if_sym_t;
  vector[N_days_tot]  rr_diag_sym_vs_sev_t;

// delay probability of reporting w x days delay (PDF of delay distribution)
  vector[Max_delay]  sym_diag_delay;
  vector[Max_delay]  sev_diag_delay;
  
  vector[Max_delay]  cas_rep_delay;
  vector[Max_delay]  die_rep_delay;
  
  vector[Max_delay]  cas_cum_report_delay; 
  vector[Max_delay]  die_cum_report_delay; 
  
// OUTCOMES
  vector[N_days_tot]  new_sym; // new in state per day
  vector[N_days_tot]  new_sev;
  vector[N_days_tot]  new_die;

  vector[N_days_tot]  new_sym_dx; 
  vector[N_days_tot]  dx_sym_sev; 
  vector[N_days_tot]  dx_sym_die; 
  
  vector[N_days_tot]  new_sev_dx;
  vector[N_days_tot]  dx_sev_die; 
  
  vector[N_days_tot]  new_die_dx;
  vector[N_days_tot]  diag_all;

  vector[N_days_tot]  occur_cas; // reported cases by date occurence
  vector[N_days_tot]  occur_die; 
  
  real                phi_cas;
  real                phi_die;

///  p_die_if_sym & p_diag_if_sym
  rr_diag_sym_vs_sev_t = inv_logit( spl_basis_dx * spl_par_dx );

  p_die_if_inf = p_sym_if_inf * p_sev_if_sym * p_die_if_sev;
  p_diag_if_sym_t = p_diag_if_sev * rr_diag_sym_vs_sev_t;

// DELAYS /////////////////////////
  for(i in 1:Max_delay){
    sym_diag_delay[i] = gamma_cdf(i+0.0, sym_prg_delay_shap, sym_prg_delay_rate/scale_dx_delay_sym)
      - gamma_cdf(i-1.0, sym_prg_delay_shap, sym_prg_delay_rate/scale_dx_delay_sym);
    sev_diag_delay[i] = gamma_cdf(i+0.0, sev_prg_delay_shap, sev_prg_delay_rate/scale_dx_delay_sev)
      - gamma_cdf(i-1.0, sev_prg_delay_shap, sev_prg_delay_rate/scale_dx_delay_sev);
  }
  
// reporting delays//~~
  for(i in 1:Max_delay) {
    cas_rep_delay[i] = gamma_cdf(i+0.0, cas_rep_delay_shap,cas_rep_delay_rate)
      - gamma_cdf(i-1.0, cas_rep_delay_shap, cas_rep_delay_rate);
    die_rep_delay[i] = gamma_cdf(i+0.0, die_rep_delay_shap, die_rep_delay_rate)
      - gamma_cdf(i-1.0, die_rep_delay_shap, die_rep_delay_rate);
}

// cumlative delays
  cas_cum_report_delay = cumulative_sum(cas_rep_delay);
  die_cum_report_delay = cumulative_sum(die_rep_delay);

// CASCADE OF INCIDENT OUTCOMES (TOTAL) ///////////////////////////////////
  new_sym = rep_vector(0, N_days_tot);
  new_sev = rep_vector(0, N_days_tot);
  new_die = rep_vector(0, N_days_tot);
  
  // 1st and 2nd difference of spl_par
  for(i in 1:(N_spl_par-2)) {
    deriv2_spl_par[i] = spl_par[i+1] * 2 - spl_par[i] - spl_par[i+2];
  }
  for(i in 1:(N_spl_par-1)) {
    deriv1_spl_par[i] = spl_par[i+1] - spl_par[i];
  }
  
  // Rt, INCIDENCE, AND NATURAL HISTORY

  logRt = spl_basis * spl_par;
  Rt  = exp(logRt);
  deriv1_log_new_inf = logRt/serial_i; 
  
  log_new_inf = cumulative_sum(deriv1_log_new_inf);
  log_new_inf = log_new_inf - log_new_inf[N_days_before] + log_new_inf_0;
  new_inf = exp(log_new_inf) + inf_imported;
  
  for(i in 1:N_days_tot) {
    for(j in 1:Max_delay) {
      if(i+(j-1) <= N_days_tot){
        new_sym[i+(j-1)] += new_inf[i] * p_sym_if_inf * inf_prg_delay[j];
      }
    }
  }

  for(i in 1:N_days_tot) {
    for(j in 1:Max_delay) {
      if(i+(j-1) <= N_days_tot){
        new_sev[i+(j-1)] += new_sym[i] * p_sev_if_sym * sym_prg_delay[j];
      }
    }
  }
  
  for(i in 1:N_days_tot) {
    for(j in 1:Max_delay) {
      if(i+(j-1) <= N_days_tot){
        new_die[i+(j-1)] += new_sev[i] * p_die_if_sev * sev_prg_delay[j];
      }
    }
  } 

// CASCADE OF INCIDENT OUTCOMES << DIAGNOSED >> ///////////
  new_sym_dx = rep_vector(0, N_days_tot);
  dx_sym_sev = rep_vector(0, N_days_tot);
  dx_sym_die = rep_vector(0, N_days_tot);
  
  new_sev_dx = rep_vector(0, N_days_tot);
  dx_sev_die = rep_vector(0, N_days_tot);
  
  diag_all   = rep_vector(0, N_days_tot);
  new_die_dx = rep_vector(0, N_days_tot);
  
// diagnosed at symptomatic
  for(i in 1:N_days_tot){
    for(j in 1:Max_delay){
      if(i+(j-1) <= N_days_tot){
        new_sym_dx[i+(j-1)] += new_sym[i] * p_diag_if_sym_t[i] * sym_diag_delay[j];
      }
    }
  }
  
// cascade from diagnosis 
  for(i in 1:N_days_tot){
    for(j in 1:Max_delay){
      if(i+(j-1) <= N_days_tot){
        dx_sym_sev[i+(j-1)] += new_sym[i] * p_diag_if_sym_t[i] * p_sev_if_sym * sym_prg_delay[j];
      }
    }
  }
  for(i in 1:N_days_tot){
    for(j in 1:Max_delay){
      if(i+(j-1) <= N_days_tot){
        dx_sym_die[i+(j-1)] += dx_sym_sev[i] * p_die_if_sev * sev_prg_delay[j];
      }
    }   
  }
    
  // diagnosed at severe 
  for(i in 1:N_days_tot){
    for(j in 1:Max_delay){
      if(i+(j-1) <= N_days_tot){
        new_sev_dx[i+(j-1)] += (new_sev[i] - dx_sym_sev[i]) * p_diag_if_sev * sev_diag_delay[j]; 
      }
    }
  }
  
  for(i in 1:N_days_tot){
    for(j in 1:Max_delay){
      if(i+(j-1) <= N_days_tot){
        dx_sev_die[i+(j-1)] += (new_sev[i] - dx_sym_sev[i]) * p_diag_if_sev * p_die_if_sev * sev_prg_delay[j];
      }
    }
  }

 //DIAGNOSIS //
  diag_all = new_sym_dx + new_sev_dx;
  new_die_dx = dx_sym_die + dx_sev_die;

// phi
  phi_cas = pow(inv_sqrt_phi_c,-2);
  phi_die = pow(inv_sqrt_phi_d,-2);
  
// REPORTING //~~
  occur_cas = rep_vector(0, N_days_tot);
  occur_die = rep_vector(0, N_days_tot);

  if(obs_cas_rep == 1) {
    for(i in 1:N_days_tot) {
      for(j in 1:Max_delay) {
        if(i+(j-1) <= N_days_tot) {
          occur_cas[i+(j-1)] += diag_all[i] * cas_rep_delay[j];
        }
      }
    }  
  } else {
    for(i in 1:(N_days_tot - Max_delay )){
      occur_cas[i] = diag_all[i]; 
    }
    for(i in (N_days_tot - Max_delay + 1):N_days_tot)  {
      occur_cas[i] += diag_all[i] * cas_cum_report_delay[N_days_tot - i + 1];
    }
  }

  if(obs_die_rep == 1) {
    for(i in 1:N_days_tot){
      for(j in 1:Max_delay){
        if(i+(j-1) <= N_days_tot){
          occur_die[i+(j-1)] += new_die_dx[i] * die_rep_delay[j];
        }
      }
    }
  } else {
    for(i in 1:(N_days_tot - Max_delay)){
      occur_die[i] = new_die_dx[i]; 
    }
    for(i in (N_days_tot - Max_delay + 1):N_days_tot)  {
      occur_die[i] += new_die_dx[i] * cas_cum_report_delay[N_days_tot - i + 1];
    }
  }

}
///////////////////////////////////////////////////////////  
model {
  int  tmp_obs_cas;
  real tmp_occur_cas;
  int  tmp_obs_die;
  real tmp_occur_die;
//// PRIORS
  // INCIDENCE
  inf_imported        ~ normal(pri_inf_imported_mu,pri_inf_imported_sd);
  log_new_inf_0       ~ normal(pri_log_new_inf_0_mu, pri_log_new_inf_0_sd);
  spl_par             ~ normal(pri_logRt_mu, pri_logRt_sd);
  serial_i            ~ lognormal(pri_serial_i_a, pri_serial_i_b);
  deriv1_spl_par      ~ student_t(10, 0, pri_deriv1_spl_par_sd);
  deriv2_spl_par      ~ student_t(10, 0, pri_deriv2_spl_par_sd);

  // SYMPTOMS AND CARE
  p_sym_if_inf        ~ beta(pri_p_sym_if_inf_a, pri_p_sym_if_inf_b);
  p_sev_if_sym        ~ beta(pri_p_sev_if_sym_a, pri_p_sev_if_sym_b);
  p_die_if_sev        ~ beta(pri_p_die_if_sev_a, pri_p_die_if_sev_b);
  p_die_if_inf        ~ beta(pri_p_die_if_inf_a, pri_p_die_if_inf_b);

  scale_dx_delay_sym  ~ beta(scale_dx_delay_sym_a, scale_dx_delay_sym_b); 
  scale_dx_delay_sev  ~ beta(scale_dx_delay_sev_a, scale_dx_delay_sev_b);
  
  cas_rep_delay_shap  ~ lognormal(pri_cas_rep_delay_shap_a, pri_cas_rep_delay_shap_b);  
  cas_rep_delay_rate  ~ lognormal(pri_cas_rep_delay_rate_a, pri_cas_rep_delay_rate_b); 
  die_rep_delay_shap  ~ lognormal(pri_die_rep_delay_shap_a, pri_die_rep_delay_shap_b);  
  die_rep_delay_rate  ~ lognormal(pri_die_rep_delay_rate_a, pri_die_rep_delay_rate_b);

// DIAGNOSIS    
  inv_logit(spl_par_dx) ~ beta(pri_rr_diag_sym_vs_sev_a, pri_rr_diag_sym_vs_sev_b);
  p_diag_if_sev       ~ beta(pri_p_diag_if_sev_a, pri_p_diag_if_sev_b);
  inv_sqrt_phi_c      ~ normal(0, pri_inv_sqrt_phi);
  inv_sqrt_phi_d      ~ normal(0, pri_inv_sqrt_phi);
    
//// LIKELIHOOD
  if(cas_yes==1){
    tmp_obs_cas = obs_cas[1];
    tmp_occur_cas = occur_cas[1 + N_days_before];
    for(i in 1:N_days) {
      target += neg_binomial_2_lpmf(tmp_obs_cas | tmp_occur_cas, phi_cas)/n_day_av;
      if(i>nda0){
        tmp_obs_cas   -= obs_cas[i - nda0];
        tmp_occur_cas -= occur_cas[i + N_days_before - nda0];
      }
      if(i<N_days){
        tmp_obs_cas   += obs_cas[i + 1];
        tmp_occur_cas += occur_cas[i + N_days_before + 1];
      }
    }
  }
  if(die_yes==1){
    tmp_obs_die = obs_die[1];
    tmp_occur_die = occur_die[1 + N_days_before];
    for(i in 1:N_days) {
      target += neg_binomial_2_lpmf(tmp_obs_die | tmp_occur_die, phi_die)/n_day_av;
      if(i>nda0){  
        tmp_obs_die   -= obs_die[i - nda0];
        tmp_occur_die -= occur_die[i + N_days_before - nda0];
      }
      if(i<N_days){
        tmp_obs_die   += obs_die[i + 1];
        tmp_occur_die += occur_die[i + N_days_before + 1];
      }
    }
  }
}
///////////////////////////////////////////////////////////
generated quantities {

}

///////////////////////////////////////////
transformed parameters {
///~~~~~~~ Define ~~~~~~~
// INCIDENCE
  vector[N_days_tot]      log_new_inf;
  vector[N_days_tot]      new_inf;

  vector[N_days_tot-2]    deriv2_log_new_inf;
  
// DELAYS probabilitily of reporting w x days delay (PDF of delay distribution)
  real<lower=0>  inf_prg_delay_rate; //~~
  real<lower=0>  sym_prg_delay_rate; //~~
  real<lower=0>  hos_prg_delay_rate; //~~
  
  real<lower=0>  inf_res_delay_rate; //~~
  real<lower=0>  sym_res_delay_rate; //~~
  real<lower=0>  hos_res_delay_rate; //~~
  
  real<lower=0>  cas_rep_delay_rate; //~~
  real<lower=0>  hos_rep_delay_rate; //~~
  real<lower=0>  die_rep_delay_rate; //~~

  vector[N_days_tot]  inf_prg_delay;
  vector[N_days_tot]  sym_prg_delay;
  vector[N_days_tot]  hos_prg_delay;
  
  vector[N_days_tot]  inf_res_delay;
  vector[N_days_tot]  sym_res_delay;
  vector[N_days_tot]  hos_res_delay;
  
  vector[N_days_tot]  cas_rep_delay;
  vector[N_days_tot]  hos_rep_delay;
  vector[N_days_tot]  die_rep_delay;
  
  // cumulative report delays
  vector[N_days_tot]  cas_cum_report_delay; 
  vector[N_days_tot]  hos_cum_report_delay; 
  vector[N_days_tot]  die_cum_report_delay; 

// OUTCOMES
  vector[N_days_tot]  new_sym; // new in state per day
  vector[N_days_tot]  new_hos;
  vector[N_days_tot]  new_die;
  
  vector[N_days_tot]  res_inf; // resolve from state (all back to no disease)
  vector[N_days_tot]  res_sym;
  vector[N_days_tot]  res_hos;

  vector[N_days_tot]  new_sym_u; // enter new state << undiagnosed >>
  vector[N_days_tot]  new_hos_u;
  vector[N_days_tot]  new_die_u;
  
  vector[N_days_tot]  res_inf_u; // resolve from state << undiagnosed >>
  vector[N_days_tot]  res_sym_u;
  vector[N_days_tot]  res_hos_u;
  
  vector[N_days_tot]  cur_inf_u; // currently in state << undiagnosed >>
  vector[N_days_tot]  cur_sym_u;
  vector[N_days_tot]  cur_hos_u;

  vector[N_days_tot]  diag_inf; // incident diagnosis
  vector[N_days_tot]  diag_sym;
  vector[N_days_tot]  diag_hos;
  vector[N_days_tot]  diag_all;
  
  vector[N_days_tot]  occur_cas; // reported cases by date occurence
  vector[N_days_tot]  occur_hos; 
  vector[N_days_tot]  occur_die; 
  
  vector[N_days_tot]  repor_cas; // reported cases by date report
  vector[N_days_tot]  repor_hos; 
  vector[N_days_tot]  repor_die; 
  
  real                phi_cas;
  real                phi_hos;
  real                phi_die;

  phi_cas = pow(inv_sqrt_phi_c,-2);
  phi_hos = pow(inv_sqrt_phi_h,-2);
  phi_die = pow(inv_sqrt_phi_d,-2);
  
// INCIDENCE

   log_new_inf[1] = log_new_inf_0;
   for(i in 1:(N_days_tot-1)) {
    log_new_inf[i+1] =  log_new_inf[i] + deriv1_log_new_inf[i] ;
  }
  new_inf = exp(log_new_inf);
  for(i in 1:(N_days_tot-2)) {
    deriv2_log_new_inf[i] = log_new_inf[i+1] * 2 - log_new_inf[i] - log_new_inf[i+2];
  }

// DLEAYS /////////////////////////
  inf_prg_delay_rate = inf_prg_delay_shap/inf_prg_delay_mn;
  sym_prg_delay_rate = sym_prg_delay_shap/sym_prg_delay_mn;
  hos_prg_delay_rate = hos_prg_delay_shap/hos_prg_delay_mn;
  
  inf_res_delay_rate = inf_res_delay_shap/inf_res_delay_mn;
  sym_res_delay_rate = sym_res_delay_shap/sym_res_delay_mn;
  hos_res_delay_rate = hos_res_delay_shap/hos_res_delay_mn;

  cas_rep_delay_rate = cas_rep_delay_shap/cas_rep_delay_mn;
  hos_rep_delay_rate = hos_rep_delay_shap/hos_rep_delay_mn;
  die_rep_delay_rate = die_rep_delay_shap/die_rep_delay_mn;
  
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
  
// DELAYS // reporting  //~~
  for(i in 1:N_days_tot) {
    cas_rep_delay[i] = gamma_cdf(i+0.0, cas_rep_delay_shap, cas_rep_delay_rate)
      - gamma_cdf(i-1.0, cas_rep_delay_shap, cas_rep_delay_rate);
    hos_rep_delay[i] = gamma_cdf(i+0.0, hos_rep_delay_shap, hos_rep_delay_rate)
      - gamma_cdf(i-1.0, hos_rep_delay_shap, hos_rep_delay_rate);
    die_rep_delay[i] = gamma_cdf(i+0.0, die_rep_delay_shap, die_rep_delay_rate)
      - gamma_cdf(i-1.0, die_rep_delay_shap, die_rep_delay_rate);
}

cas_cum_report_delay = cumulative_sum(cas_rep_delay);
hos_cum_report_delay = cumulative_sum(hos_rep_delay);
die_cum_report_delay = cumulative_sum(die_rep_delay);
  
// CASCADE OF INCIDENT OUTCOMES (TOTAL) ///////////////////////////////////
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
        (1-p_sym_if_inf) * inf_res_delay[j]; // probability resolved on day j
      }
    }
  }
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

// CASCADE OF INCIDENT OUTCOMES << UNDIAGNOSED >> ///////////
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

// REPORTING EXPERIMENT   //~~
 // check way reporting delay is counted -- number reported on day i  
 // reporting triangle more flexible, can be done either way
 
  occur_cas = rep_vector(0, N_days_tot);
  occur_hos = rep_vector(0, N_days_tot);
  occur_die = rep_vector(0, N_days_tot);
  
// for data by diagnosis date
for(i in 1:N_days_tot)  {
  occur_cas[i] += diag_all[i] * cas_cum_report_delay[N_days_tot - i + 1];
}
  
for(i in 1:N_days_tot)  {
  occur_hos[i] += diag_hos[i] * hos_cum_report_delay[N_days_tot - i + 1];
}

for(i in 1:N_days_tot)  {
  occur_die[i] += (new_die[i] - new_die_u[i]) * die_cum_report_delay[N_days_tot - i + 1];
}

  repor_cas = rep_vector(0, N_days_tot);
  repor_hos = rep_vector(0, N_days_tot);
  repor_die = rep_vector(0, N_days_tot);

// for data by reporting date  
 for(i in 1:N_days_tot){
    repor_cas[i] += diag_all[i] * cas_rep_delay[i];
 }
 
  for(i in 1:N_days_tot){
    repor_hos[i] += diag_hos[i] * hos_rep_delay[i];
 }

 for(i in 1:N_days_tot){
    repor_die[i] += (new_die[i] - new_die_u[i]) *  die_rep_delay[i]
 }

}

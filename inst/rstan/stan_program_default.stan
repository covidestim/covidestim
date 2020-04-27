data {
  ///~~~~~~~ Define ~~~~~~~
  // INPUT DATA
  int<lower=0>     N_days;
  int<lower=0>     N_days_delay; //~~ // days BEFORE data to initialized epi model
  int<lower=0>     obs_cas[N_days]; //~~
  int<lower=0>     obs_hos[N_days]; //~~
  int<lower=0>     obs_die[N_days]; //~~

  // PRIORS

  // random walk 
  real              pri_log_new_inf_0_mu;
  real<lower=0>     pri_log_new_inf_0_sd;
  real<lower=0>     pri_deriv1_log_new_inf_sd; //~~ new name
  real<lower=0>     pri_deriv2_log_new_inf_sd;
  // p(progression)
  real<lower=0>          pri_p_sym_if_inf_a; 
  real<lower=0>          pri_p_sym_if_inf_b;
  real<lower=0>          pri_p_hos_if_sym_a;
  real<lower=0>          pri_p_hos_if_sym_b;
  real<lower=0>          pri_p_die_if_hos_a;
  real<lower=0>          pri_p_die_if_hos_b;
  // p(diag)  
  real<lower=0>          pri_p_diag_if_inf_a;
  real<lower=0>          pri_p_diag_if_inf_b;
  real<lower=0>          pri_p_diag_if_sym_a;
  real<lower=0>          pri_p_diag_if_sym_b;
  real<lower=0>          pri_p_diag_if_hos_a;
  real<lower=0>          pri_p_diag_if_hos_b;
  // delay to progression
  real<lower=0>          pri_inf_prg_delay_shap; //~~ new names
  real<lower=0>          pri_inf_prg_delay_rate;
  real<lower=0>          pri_sym_prg_delay_shap;
  real<lower=0>          pri_sym_prg_delay_rate;
  real<lower=0>          pri_hos_prg_delay_shap;
  real<lower=0>          pri_hos_prg_delay_rate;
  // delay to recovered  
  real<lower=0>          pri_inf_res_delay_shap; //~~ new names
  real<lower=0>          pri_inf_res_delay_rate;
  real<lower=0>          pri_sym_res_delay_shap;
  real<lower=0>          pri_sym_res_delay_rate;
  real<lower=0>          pri_hos_res_delay_shap;
  real<lower=0>          pri_hos_res_delay_rate;
  // delay to report
  real<lower=0>          pri_cas_rep_delay_shap;
  real<lower=0>          pri_cas_rep_delay_rate;
  real<lower=0>          pri_hos_rep_delay_shap;
  real<lower=0>          pri_hos_rep_delay_rate;
  real<lower=0>          pri_die_rep_delay_shap;
  real<lower=0>          pri_die_rep_delay_rate;
  // additional delay uncertainties
  real<lower=0>          inf_prg_delay_shap_a; 
  real<lower=0>          inf_prg_delay_shap_b;
  real<lower=0>          sym_prg_delay_shap_a;
  real<lower=0>          sym_prg_delay_shap_b;
  real<lower=0>          hos_prg_delay_shap_a;
  real<lower=0>          hos_prg_delay_shap_b;

  real<lower=0>          inf_res_delay_shap_a; 
  real<lower=0>          inf_res_delay_shap_b;
  real<lower=0>          sym_res_delay_shap_a;
  real<lower=0>          sym_res_delay_shap_b;
  real<lower=0>          hos_res_delay_shap_a;
  real<lower=0>          hos_res_delay_shap_b;

  real<lower=0>          cas_rep_delay_shp_a; 
  real<lower=0>          cas_rep_delay_shp_b;
  real<lower=0>          hos_rep_delay_shp_a; 
  real<lower=0>          hos_rep_delay_shp_b;
  real<lower=0>          die_rep_delay_shp_a; 
  real<lower=0>          die_rep_delay_shp_b;

  // turn on/off negative binomial
  int<lower = 0, upper = 1> nb_yes; 
  // set whether data are by diagnosis date (vs. reporting date)
  int<lower=0, upper=1> obs_cas_rep; //~~
  int<lower=0, upper=1> obs_hos_rep; //~~
  int<lower=0, upper=1> obs_die_rep; //~~
}
///////////////////////////////////////////////////////////
transformed data {

  int<lower=1>  N_days_tot;
  
  N_days_tot = N_days + N_days_delay; 


  // KEEP COMMENTED OUT: important for later version that take linelist data. 

  // reporting triangle -- laying out cases by date of report and by delay.
  //int           rep_tri_conf_cases[N_days, Max_delay+1]; 
  
  //N_days_tot = N_days + N_days_extra; // how many days to project for
  //for(i in 1:N_days){
  //  for(j in 1:(Max_delay+1)){
  //    rep_tri_conf_cases[i,j] = 0;
  //  }
  //}
  //for(i in 1:N_conf_cases){
  //  rep_tri_conf_cases[cases_test_day[i], cases_days_delay[i]+1] += 1;
  //}
}
///////////////////////////////////////////////////////////
parameters {
  
// INCIDENCE  
 real                  log_new_inf_0; // intercept in log space
 vector[N_days_tot-1]  deriv1_log_new_inf; // first derivative of the random walk 

// SYMPTOMS AND CARE 
  real<lower=0, upper=1>  p_sym_if_inf;
  real<lower=0, upper=1>  p_hos_if_sym;
  real<lower=0, upper=1>  p_die_if_hos;
  
// delay associated with transition to next illness state
  real<lower=1>           inf_prg_delay_mn;
  real<lower=1>           sym_prg_delay_mn;
  real<lower=1>           hos_prg_delay_mn;
// delay associated with transition to recovered 
  real<lower=1>           inf_res_delay_mn;
  real<lower=1>           sym_res_delay_mn;
  real<lower=1>           hos_res_delay_mn;
// delay associated with reporting   
  real<lower=1>           cas_rep_delay_mn; //~~         
  real<lower=1>           hos_rep_delay_mn; //~~       
  real<lower=1>           die_rep_delay_mn; //~~           
// uncertainty in delays
  real<lower=1>           inf_prg_delay_shap; //~~
  real<lower=1>           sym_prg_delay_shap; //~~
  real<lower=1>           hos_prg_delay_shap; //~~
  
  real<lower=1>           inf_res_delay_shap; //~~
  real<lower=1>           sym_res_delay_shap; //~~
  real<lower=1>           hos_res_delay_shap; //~~
  
  real<lower=1>           cas_rep_delay_shap; //~~
  real<lower=1>           hos_rep_delay_shap; //~~
  real<lower=1>           die_rep_delay_shap; //~~
    
// DIAGNOSIS // probability of diagnosis at each illness state
  real<lower=0, upper=1>  p_diag_if_inf;
  real<lower=0, upper=1>  p_diag_if_sym;
  real<lower=0, upper=1>  p_diag_if_hos;
  
// LIKELIHOOD
  real<lower=0>           inv_sqrt_phi_c;
  real<lower=0>           inv_sqrt_phi_h;
  real<lower=0>           inv_sqrt_phi_d;

}

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
  
  vector[N_days_tot]  cur_inf; // currently in state << undiagnosed >>
  vector[N_days_tot]  cur_sym;
  vector[N_days_tot]  cur_hos;

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
    for(j in 1:(N_days_tot-i+1)) {  
      res_inf[i+(j-1)] += new_inf[i] * (1-p_sym_if_inf) *  inf_res_delay[j]; 
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

// REPORTING //~~

  occur_cas = rep_vector(0, N_days_tot);
  occur_hos = rep_vector(0, N_days_tot);
  occur_die = rep_vector(0, N_days_tot);
  
// for data by diagnosis date
for(i in 1:N_days_tot)  {
  occur_cas[i] += diag_all[i] * cas_cum_report_delay[N_days_tot - i + 1];
  occur_hos[i] += (new_hos[i] - new_hos_u[i] + diag_hos[i]) * hos_cum_report_delay[N_days_tot - i + 1];
  occur_die[i] += (new_die[i] - new_die_u[i]) * die_cum_report_delay[N_days_tot - i + 1];
}

  repor_cas = rep_vector(0, N_days_tot);
  repor_hos = rep_vector(0, N_days_tot);
  repor_die = rep_vector(0, N_days_tot);

// for data by reporting date  
for(i in 1:N_days_tot){
   for(j in 1:(N_days_tot - i +1)){
    repor_cas[i+(j-1)] += diag_all[i] * cas_rep_delay[j];
  }
}
  for(i in 1:N_days_tot){
    for(j in 1:(N_days_tot - i +1)){
    repor_hos[i+(j-1)] += (cur_sym[i] - cur_sym_u[i] + diag_hos[i]) * 
      hos_rep_delay[j];
  }
}

 for(i in 1:N_days_tot){
   for(j in 1:(N_days_tot - i +1)){
    repor_die[i+(j-1)] += (new_die[i] - new_die_u[i]) *  die_rep_delay[j];
  }
}  

}
///////////////////////////////////////////////////////////  
model {
//// PRIORS
  // INCIDENCE
   log_new_inf_0             ~ normal(pri_log_new_inf_0_mu, pri_log_new_inf_0_sd); // 
   deriv1_log_new_inf        ~ normal(0, pri_deriv1_log_new_inf_sd); //~~
   deriv2_log_new_inf        ~ normal(0, pri_deriv2_log_new_inf_sd);
  
  // SYMPTOMS AND CARE
    p_sym_if_inf              ~ beta(pri_p_sym_if_inf_a, pri_p_sym_if_inf_b);
    p_hos_if_sym              ~ beta(pri_p_hos_if_sym_a, pri_p_hos_if_sym_b);
    p_die_if_hos              ~ beta(pri_p_die_if_hos_a, pri_p_die_if_hos_b);
    
    inf_prg_delay_shap        ~ normal(inf_prg_delay_shap_a, inf_prg_delay_shap_b); //~~
    sym_prg_delay_shap        ~ normal(sym_prg_delay_shap_a, sym_prg_delay_shap_b); //~~
    hos_prg_delay_shap        ~ normal(hos_prg_delay_shap_a, hos_prg_delay_shap_b); //~~
    inf_res_delay_shap        ~ normal(inf_res_delay_shap_a, inf_res_delay_shap_b); //~~
    sym_res_delay_shap        ~ normal(sym_res_delay_shap_a, sym_res_delay_shap_b); //~~
    hos_res_delay_shap        ~ normal(hos_res_delay_shap_a, hos_res_delay_shap_b); //~~
    
    cas_rep_delay_shap        ~ normal(cas_rep_delay_shp_a, cas_rep_delay_shp_b); //~~
    hos_rep_delay_shap        ~ normal(hos_rep_delay_shp_a, hos_rep_delay_shp_b); //~~
    die_rep_delay_shap        ~ normal(die_rep_delay_shp_a, die_rep_delay_shp_b); //~~
    
    // new priors names
    inf_prg_delay_mn          ~ gamma(pri_inf_prg_delay_shap, pri_inf_prg_delay_rate); 
    sym_prg_delay_mn          ~ gamma(pri_sym_prg_delay_shap, pri_sym_prg_delay_rate);
    hos_prg_delay_mn          ~ gamma(pri_hos_prg_delay_shap, pri_hos_prg_delay_rate);
    inf_res_delay_mn          ~ gamma(pri_inf_res_delay_shap, pri_inf_res_delay_rate);
    sym_res_delay_mn          ~ gamma(pri_sym_res_delay_shap, pri_sym_res_delay_rate);
    hos_res_delay_mn          ~ gamma(pri_hos_res_delay_shap, pri_hos_res_delay_rate);
    
    // add priors for each delay dist
    cas_rep_delay_mn          ~gamma(pri_cas_rep_delay_shap, pri_cas_rep_delay_rate); //~~
    hos_rep_delay_mn          ~gamma(pri_hos_rep_delay_shap, pri_hos_rep_delay_rate); //~~
    die_rep_delay_mn          ~gamma(pri_die_rep_delay_shap, pri_die_rep_delay_rate); //~~
  
    p_diag_if_inf             ~ beta(pri_p_diag_if_inf_a, pri_p_diag_if_inf_b);
    p_diag_if_sym             ~ beta(pri_p_diag_if_sym_a, pri_p_diag_if_sym_b);
    p_diag_if_hos             ~ beta(pri_p_diag_if_hos_a, pri_p_diag_if_hos_b);
    
    inv_sqrt_phi_c            ~ normal(0, 1); //~~
    inv_sqrt_phi_h            ~ normal(0, 1); //~~
    inv_sqrt_phi_d            ~ normal(0, 1); //~~
    
//// LIKELIHOOD
//SWITCH TO NEG BIN
  if (nb_yes == 1) {
      if(obs_cas_rep == 0){
         for(i in 1:N_days) {
            obs_cas[i] ~ neg_binomial_2(occur_cas[i + N_days_delay], phi_cas);
          }
      } else {
         for(i in 1:N_days) {
            obs_cas[i] ~ neg_binomial_2(repor_cas[i + N_days_delay], phi_cas);
          }
       }
      
      if(obs_hos_rep == 0){
         for(i in 1:N_days) {
            obs_hos[i] ~ neg_binomial_2(occur_hos[i + N_days_delay], phi_hos);
              }
          } else {
         for(i in 1:N_days) {
            obs_hos[i] ~ neg_binomial_2(repor_hos[i + N_days_delay], phi_hos);
           }
        }
     
     if(obs_die_rep == 0){
       for(i in 1:N_days) {
          obs_die[i] ~ neg_binomial_2(occur_die[i + N_days_delay], phi_die);
         }
      } else {
      for(i in 1:N_days) {
          obs_die[i] ~ neg_binomial_2(repor_die[i + N_days_delay], phi_die);
        }
     }
     
  } else { 
       if(obs_cas_rep == 0){
         for(i in 1:N_days) {
            obs_cas[i] ~ poisson(occur_cas[i + N_days_delay]);
          }
      } else {
         for(i in 1:N_days) {
            obs_cas[i] ~ poisson(repor_cas[i + N_days_delay]);
          }
      }
      if(obs_hos_rep == 0){
         for(i in 1:N_days) {
            obs_hos[i] ~ poisson(occur_hos[i + N_days_delay]);
            }
       } else {
         for(i in 1:N_days) {
            obs_hos[i] ~ poisson(repor_hos[i + N_days_delay]);
          }
      }
     if(obs_die_rep == 0){
       for(i in 1:N_days) {
          obs_die[i] ~ poisson(occur_die[i + N_days_delay]);
          }
      } else {
      for(i in 1:N_days) {
          obs_die[i] ~ poisson(repor_die[i + N_days_delay]);
        }
     }
   }
}
///////////////////////////////////////////////////////////
generated quantities {
}

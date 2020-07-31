data {
  ///~~~~~~~ Define ~~~~~~~
  // INPUT DATA
  int<lower=0>           N_days; // days of data
  int<lower=0>           N_days_before; // days before data to initialized epi model
  int<lower=0>           Max_delay; // maximum days delay 
  
  int<lower=0>           obs_cas[N_days]; // vector of cases
  int<lower=0>           obs_die[N_days]; // vector of deaths
  real<lower=0,upper=1>  frac_pos[N_days+N_days_before]; // vector of fraction 
                                                         //positive tests
  // modulate relationship between testing fraction, diagnosis
  real<lower=0,upper=1>  rho_sym; 
  real<lower=0,upper=1>  rho_sev; 
  int<lower=0,upper=1>   is_weekend[N_days+N_days_before]; // weekend indicator
  
  // terms for spline
  int<lower=1> n_spl_par;
  matrix[(N_days + N_days_before), n_spl_par] spl_basis;
  
  // delay distributions! 
  // fixed delay distribtions. time from inf -> sym, sym -> sev, sev -> die
  real<lower=0>          inf_prg_delay_shap; 
  real<lower=0>          inf_prg_delay_rate;
  real<lower=0>          sym_prg_delay_shap;
  real<lower=0>          sym_prg_delay_rate;
  real<lower=0>          sev_prg_delay_shap;
  real<lower=0>          sev_prg_delay_rate;
  // fixed delay from diagnosis to report
  real<lower=0>          cas_rep_delay_shap;
  real<lower=0>          cas_rep_delay_rate;
  real<lower=0>          die_rep_delay_shap;
  real<lower=0>          die_rep_delay_rate;
  
  /// control nobs! 
  //  what data are included --  cases, deaths: 
  int<lower = 0, upper = 1> cas_yes; 
  int<lower = 0, upper = 1> die_yes; 
  //  how are the data dated -- report, occurrence: 
  int<lower = 0, upper = 1> obs_cas_rep;
  int<lower = 0, upper = 1> obs_die_rep;
  //  how many days should be used for the moving average in the likelihood 
  //  function? 
  int<lower = 1, upper = 10> N_days_av; 

  // TERMS FOR PRIOR DISTRIBTUIONS
  // for new infections
  real                   pri_log_new_inf_0_mu;
  real<lower=0>          pri_log_new_inf_0_sd;
  real                   pri_logRt_mu; // ~~ NEW
  real<lower=0>          pri_logRt_sd; // ~~ NEW
  real                   pri_serial_i_a; // ~~ NEW
  real<lower=0>          pri_serial_i_b; // ~~ NEW
  real                   pri_inf_imported_mu; // ~~ NEW
  real<lower=0>          pri_inf_imported_sd; // ~~ NEW
  real<lower=0>          pri_deriv1_b_spline_sd; // ~~ NEW
  real<lower=0>          pri_deriv2_b_spline_sd; // ~~ NEW
  // probabilities of progression inf -> sym -> sev -> die
  real<lower=0>          pri_p_sym_if_inf_a; 
  real<lower=0>          pri_p_sym_if_inf_b;
  real<lower=0>          pri_p_sev_if_sym_a;
  real<lower=0>          pri_p_sev_if_sym_b;
  real<lower=0>          pri_p_die_if_sev_a;
  real<lower=0>          pri_p_die_if_sev_b;
  // overall case fatality rate
  real<lower=0>          pri_p_die_if_sym_a;
  real<lower=0>          pri_p_die_if_sym_b;
  // probabilities of diagnosis 
  real<lower=0>          pri_p_diag_if_sym_a;
  real<lower=0>          pri_p_diag_if_sym_b;
  real<lower=0>          pri_p_diag_if_sev_a;
  real<lower=0>          pri_p_diag_if_sev_b;
  // size of the 'weekend effect' on diagnosis 
  real<lower=0>          pri_weekend_eff_a; //~~    
  real<lower=0>          pri_weekend_eff_b; //~~
  // delay to diagnosis assumed to be some fraction of progression delay
  // Beta prior distribtuion for that fraction 
  real<lower=0>          scale_dx_delay_sym_a; 
  real<lower=0>          scale_dx_delay_sym_b; 
  real<lower=0>          scale_dx_delay_sev_a; 
  real<lower=0>          scale_dx_delay_sev_b;

}
///////////////////////////////////////////////////////////
transformed data {
 int  N_days_tot;
 int  nda0;
 
 // Progression delays
 vector[Max_delay]  inf_prg_delay;
 vector[Max_delay]  sym_prg_delay;
 vector[Max_delay]  sev_prg_delay;
// Reporting delays
 vector[Max_delay]  cas_rep_delay;
 vector[Max_delay]  die_rep_delay;
// Cumulative reporting delays
 vector[Max_delay]  cas_cum_report_delay; 
 vector[Max_delay]  die_cum_report_delay; 
 
// for likelihood function moving avg. 
  nda0 = N_days_av - 1; 
   
// create 'N_days_tot', which is days of data plus days to model before first 
// case or death 
  N_days_tot = N_days + N_days_before; 
  
 // calculate the daily probability of transitioning to a new disease state
 // for days 1 to 60 after entering that state
  for(i in 1:Max_delay) {
    inf_prg_delay[i] = gamma_cdf(i+0.0, inf_prg_delay_shap, inf_prg_delay_rate)
      - gamma_cdf(i-1.0, inf_prg_delay_shap, inf_prg_delay_rate);
    sym_prg_delay[i] = gamma_cdf(i+0.0, sym_prg_delay_shap, sym_prg_delay_rate)
      - gamma_cdf(i-1.0, sym_prg_delay_shap, sym_prg_delay_rate);
    sev_prg_delay[i] = gamma_cdf(i+0.0, sev_prg_delay_shap, sev_prg_delay_rate)
      - gamma_cdf(i-1.0, sev_prg_delay_shap, sev_prg_delay_rate);
  }
  
// Calcluate the probability of reporting for each day after diagnosis
// for 1 to 60 post diagnosis. 
  for(i in 1:Max_delay) {
    cas_rep_delay[i] = gamma_cdf(i+0.0, cas_rep_delay_shap,cas_rep_delay_rate)
      - gamma_cdf(i-1.0, cas_rep_delay_shap, cas_rep_delay_rate);
    die_rep_delay[i] = gamma_cdf(i+0.0, die_rep_delay_shap, die_rep_delay_rate)
      - gamma_cdf(i-1.0, die_rep_delay_shap, die_rep_delay_rate);
}

// Cumlative Reporting Delays
// Calculate the cumulative probability of reporting for each day 
// after diagnosis.
cas_cum_report_delay = cumulative_sum(cas_rep_delay);
die_cum_report_delay = cumulative_sum(die_rep_delay);
}
///////////////////////////////////////////////////////////
parameters {
  
// INCIDENCE 
 real                  log_new_inf_0; // infection intercept in log space
 vector[n_spl_par]     b_spline; // spline for Rt
 real<lower=0>         serial_i; // serial interval
 real<lower=0>         inf_imported; // imported cases

// DISEASE PROGRESSION
// probability of transitioning between disease states
  real<lower=0, upper=1>    p_sym_if_inf;
  real<lower=0, upper=1>    p_sev_if_sym;
  real<lower=0, upper=0.3>  p_die_if_sev;
  
// DIANGOSIS
// scaling factor for time to diagnosis
  real<lower=0, upper=1>    scale_dx_delay_sym; 
  real<lower=0, upper=1>    scale_dx_delay_sev; 
// probability of diagnosis at each illness state
  real<lower=0, upper=1>    p_diag_if_sym;
  real<lower=0, upper=1>    p_diag_if_sev;
  real<lower=0, upper=1>    weekend_eff;
 
// LIKELIHOOD 
// phi terms for negative b ino imal likelihood function 
  real<lower=0>             inv_sqrt_phi_c;
  real<lower=0>             inv_sqrt_phi_d;

}

///////////////////////////////////////////
transformed parameters {
///~~~~~~~ Define ~~~~~~~
// INCIDENCE
  vector[N_days_tot]      log_new_inf;
  vector[N_days_tot]      new_inf;
  vector[N_days_tot]      deriv1_log_new_inf;
  
  // Rt spline
  vector[N_days_tot]      logRt;
  vector[N_days_tot]      Rt;
  vector[n_spl_par-2]     deriv2_b_spline;
  vector[n_spl_par-1]     deriv1_b_spline;
  
// DIAGNOSIS AND REPORTING  
 // daily probabilities of diagnosis and report
 // for days 1 to 60 after entering that state
  vector[Max_delay]  sym_diag_delay;
  vector[Max_delay]  sev_diag_delay;

// DISEASE OUTCOMES
  // overall case fatality rate
  real<lower=0, upper=0.1>  p_die_if_sym;
  
  // "true" number entering disease state each day
  vector[N_days_tot]  new_sym; 
  vector[N_days_tot]  new_sev;
  vector[N_days_tot]  new_die;
  // newly diagnosed
  vector[N_days_tot]  new_sym_dx; 
  vector[N_days_tot]  new_sev_dx;
  // follow diagnosed cases forward to calculate deaths among diagnosed
  vector[N_days_tot]  dx_sym_sev; 
  vector[N_days_tot]  dx_sym_die; 
  vector[N_days_tot]  dx_sev_die; 
  // sum to diagnosed cases and deaths
  vector[N_days_tot]  diag_all;
  vector[N_days_tot]  new_die_dx;
  // number of cases and deaths in official record on each day
  // (all diagnosed cases with an additional delay to report) 
  vector[N_days_tot]  occur_cas;
  vector[N_days_tot]  occur_die; 
  
// LIKELIHOOD
// phi terms for negative bino imal likelihood function 
  real                phi_cas;
  real                phi_die;

// DLEAYS //

// Diagnosis Delays
// Calcluate the probability of diagnosis for each day in state
// for 1 to 60 days in state. We do this by scaling the rate term
// of the gamma distribution of progressiong delays by a modeled fraction 
// (scale_dx_delay_xxx)
for(i in 1:Max_delay){
    sym_diag_delay[i] = gamma_cdf(i+0.0, sym_prg_delay_shap, 
                        sym_prg_delay_rate/scale_dx_delay_sym)
                      - gamma_cdf(i-1.0, sym_prg_delay_shap, 
                      sym_prg_delay_rate/scale_dx_delay_sym);
    sev_diag_delay[i] = gamma_cdf(i+0.0, sev_prg_delay_shap, 
                        sev_prg_delay_rate/scale_dx_delay_sev)
                        - gamma_cdf(i-1.0, sev_prg_delay_shap, 
                        sev_prg_delay_rate/scale_dx_delay_sev);
}
  
// CASCADE OF INCIDENT OUTCOMES ("TRUE") //
// initialize empty vectors
  new_sym = rep_vector(0, N_days_tot);
  new_sev = rep_vector(0, N_days_tot);
  new_die = rep_vector(0, N_days_tot);
  
// NEW INCIDENT CASES
  
  // modeled with a spline
  logRt = spl_basis * b_spline;
  Rt = exp(logRt); 
  deriv1_log_new_inf = logRt/serial_i; 
  
  log_new_inf = cumulative_sum(deriv1_log_new_inf);
  log_new_inf = log_new_inf - log_new_inf[N_days_before] + log_new_inf_0;
  new_inf = exp(log_new_inf) + inf_imported;
  
  // second derivative
  for(i in 1:(n_spl_par-2)) {
    deriv2_b_spline[i] = b_spline[i+1] * 2 - b_spline[i] - b_spline[i+2];
  }
  // first derivative
  for(i in 1:(n_spl_par-1)) {
    deriv1_b_spline[i] = b_spline[i+1] - b_spline[i];
  }

 // SYMPTOMATIC CASES
 // cases entering a state on day i + j - 1: 
 // cases entering previous state on day i * the probability of progression *
 // the probability progression occurred on day j 
 
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
  
// overall case fatality rate is equal to the probability of death among
// severely ill individuals times the probability of being severely ill. 
p_die_if_sym = p_die_if_sev * p_sev_if_sym; 

// CASCADE OF INCIDENT OUTCOMES (DIAGNOSED) //
// initialize empty vectors
  new_sym_dx = rep_vector(0, N_days_tot);
    dx_sym_sev = rep_vector(0, N_days_tot);
    dx_sym_die = rep_vector(0, N_days_tot);
  
  new_sev_dx = rep_vector(0, N_days_tot);
    dx_sev_die = rep_vector(0, N_days_tot);
  
  diag_all   = rep_vector(0, N_days_tot);
  new_die_dx = rep_vector(0, N_days_tot);
  
// diagnosed at symptomatic
// a diagnosed symptomatic (not severe) case on day i + j - 1
// is a symptomatic case on day i with some probability of diagnosis
// (impacted by the fraction positive on day i and whether day i + j -1
// is a weekend), and some probability that the diagnosis occurred on day j.  
  for(i in 1:N_days_tot){
    for(j in 1:Max_delay){
      if(i+(j-1) <= N_days_tot){
      new_sym_dx[i+(j-1)] += new_sym[i]  * 
                             p_diag_if_sym * pow(1-frac_pos[i], rho_sym) * 
                             (1 - (is_weekend[i+(j-1)] * weekend_eff)) * 
                             sym_diag_delay[j];
    }
  }
}  
// cascade from diagnosis 
// follow diagnosed cases forward to determine how many cases diagnosed
// at symptomatic eventually die 
          for(i in 1:N_days_tot){
            for(j in 1:Max_delay){
              if(i+(j-1) <= N_days_tot){
               dx_sym_sev[i+(j-1)] += new_sym[i] * p_diag_if_sym * 
               pow(1-frac_pos[i], rho_sym) * 
              // (1 - (is_weekend[i+]* weekend_eff)) * 
               p_sev_if_sym * sym_prg_delay[j];
            }
          }
      }
      
      for(i in 1:N_days_tot){
            for(j in 1:Max_delay){
              if(i+(j-1) <= N_days_tot){
              dx_sym_die[i+(j-1)] += dx_sym_sev[i] * p_die_if_sev *
                sev_prg_delay[j];
            }
          }    
        }
        
// diagnosed at severe 
// as above for symptomatic 
  for(i in 1:N_days_tot){
    for(j in 1:Max_delay){
      if(i+(j-1) <= N_days_tot){
      new_sev_dx[i+(j-1)] += (new_sev[i] - dx_sym_sev[i]) * 
                             p_diag_if_sev * pow(1-frac_pos[i], rho_sev) * 
                             (1 - (is_weekend[i+(j-1)] * weekend_eff)) * 
                             sev_diag_delay[j]; 
    }
  }
} 
// cascade from diagnosis
// as above for symptomatic 
          for(i in 1:N_days_tot){
            for(j in 1:Max_delay){
              if(i+(j-1) <= N_days_tot){
              dx_sev_die[i+(j-1)] += (new_sev[i] - dx_sym_sev[i]) * 
              p_diag_if_sev *  pow(1-frac_pos[i], rho_sev) * 
              //(1 - (is_weekend[i+(j-1)] * weekend_eff)) * 
              p_die_if_sev * sev_prg_delay[j];
            }
          }
        }  

// TOTAL DIAGNOSED CASES AND DEATHS //
diag_all = new_sym_dx + new_sev_dx;
new_die_dx = dx_sym_die + dx_sev_die;

// REPORTING //
// calcluate "occur_cas" and "occur_die", which are vectors of diagnosed cases 
// and deaths by the date we expect them to appear in the reported data. 

// initialize empty vectors
occur_cas = rep_vector(0, N_days_tot);
occur_die = rep_vector(0, N_days_tot);

// how reporting delays are reflected in the data depend on how the data are 
// dated
if(obs_cas_rep == 1) {
  // for cases by date of report: 
  for(i in 1:N_days_tot) {
    for(j in 1:Max_delay) {
      if(i+(j-1) <= N_days_tot) {
  // a reported case on day i + j - 1 is a case that was diagnosed on day i and
  // has some probability of being reported on day j. 
        occur_cas[i+(j-1)] += diag_all[i] * cas_rep_delay[j];
      }
    }
  }  
} else {
  // for cases by date of occurrence
  // we assume all cases diagnosed more than 60 days from the final day of data
  // have been reported
   for(i in 1:(N_days_tot - Max_delay )){
        occur_cas[i] = diag_all[i]; 
      }
  // we model completeness of reporting with a cumulative probability of 
  // reporting on day j
  //e.g. diagnosed cases * probability of reporting 1 day after diagnosis
  //     diagnosed cases * probability of reporting 2 days after diagnosis
  // etc. such that cases become less complete as we approach the end of 
  // the time series of reported cases. 
      for(i in (N_days_tot - Max_delay + 1):N_days_tot)  {
        occur_cas[i] += diag_all[i] * cas_cum_report_delay[N_days_tot - i + 1];
      }
}
// reporting delays modeled as described above for cases
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

// phi
phi_cas = pow(inv_sqrt_phi_c, -2);
phi_die = pow(inv_sqrt_phi_d, -2);

}
///////////////////////////////////////////////////////////  
model {
  // placeholders to allow moving average in likelihood
  int  tmp_obs_cas;
  real tmp_occur_cas;
  int  tmp_obs_die;
  real tmp_occur_die;
//// PRIORS
  log_new_inf_0         ~ normal(pri_log_new_inf_0_mu, pri_log_new_inf_0_sd);
  b_spline              ~ normal(pri_logRt_mu, pri_logRt_sd);
  serial_i              ~ lognormal(pri_serial_i_a, pri_serial_i_b);
  deriv2_b_spline       ~ normal(0, pri_deriv1_b_spline_sd);
  deriv1_b_spline       ~ normal(0, pri_deriv2_b_spline_sd);
  
  // DISEASE PROGRESSION
  // prbability of transitioning from inf -> sym -> sev -> die
  p_sym_if_inf         ~ beta(pri_p_sym_if_inf_a, pri_p_sym_if_inf_b);
  p_sev_if_sym         ~ beta(pri_p_sev_if_sym_a, pri_p_sev_if_sym_b);
  p_die_if_sev         ~ beta(pri_p_die_if_sev_a, pri_p_die_if_sev_b);
  // overall case fatality rate
  p_die_if_sym         ~ beta(pri_p_die_if_sym_a, pri_p_die_if_sym_b);

// DIAGNOSIS    
  // probabilities of diagnosis
  p_diag_if_sym        ~ beta(pri_p_diag_if_sym_a, pri_p_diag_if_sym_b);
  p_diag_if_sev        ~ beta(pri_p_diag_if_sev_a, pri_p_diag_if_sev_b);
  weekend_eff          ~ beta(pri_weekend_eff_a, pri_weekend_eff_b);
  // delay distribution scaling factors
  scale_dx_delay_sym   ~ beta(scale_dx_delay_sym_a, scale_dx_delay_sym_b); 
  scale_dx_delay_sev   ~ beta(scale_dx_delay_sev_a, scale_dx_delay_sev_b);

// phi  
  inv_sqrt_phi_c       ~ normal(0, 1); //~~
  inv_sqrt_phi_d       ~ normal(0, 1); //~~
    
//// LIKELIHOOD
  if(cas_yes==1){ // if cases data are present
    tmp_obs_cas = obs_cas[1]; // cases, first day model
    tmp_occur_cas = occur_cas[1 + N_days_before]; // cases, first day data
    for(i in 1:N_days) {
      target += neg_binomial_2_lpmf(tmp_obs_cas | tmp_occur_cas, phi_cas)/N_days_av;
      if(i>nda0){ // if day is greater than or equal to moving avg. size
        tmp_obs_cas   -= obs_cas[i - nda0];
        tmp_occur_cas -= occur_cas[i + N_days_before - nda0];
      }
      if(i<N_days){ // if it is not the final day of data
        tmp_obs_cas   += obs_cas[i + 1];
        tmp_occur_cas += occur_cas[i + N_days_before + 1];
      }
    }
  }
  if(die_yes==1){ // if death data are present
    tmp_obs_die = obs_die[1];
    tmp_occur_die = occur_die[1 + N_days_before];
    for(i in 1:N_days) {
      target += neg_binomial_2_lpmf(tmp_obs_die | tmp_occur_die, phi_die)/N_days_av;
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
  // calculate cumulative incidence as the cumunlative sum of 
  // 'true' incident cases
  vector[N_days_tot]  cumulative_incidence;  
  cumulative_incidence = cumulative_sum(new_inf); 
}



data {
  ///~~~~~~~ Define ~~~~~~~
  // INPUT DATA
  int<lower=0>           N_days;
  int<lower=0>           N_days_before; //~~ // days before data to initialized epi model
  
  int<lower=0>           obs_cas[N_days]; //~~
  int<lower=0>           obs_die[N_days]; //~~
  // PRIORS
  // random walk 
  real                   pri_log_new_inf_0_mu;
  real<lower=0>          pri_log_new_inf_0_sd;
  real<lower=0>          pri_deriv1_log_new_inf_sd; //~~ new name
  real<lower=0>          pri_deriv2_log_new_inf_sd;
  // p(progression)
  real<lower=0>          pri_p_sym_if_inf_a; 
  real<lower=0>          pri_p_sym_if_inf_b;
  real<lower=0>          pri_p_sev_if_sym_a;
  real<lower=0>          pri_p_sev_if_sym_b;
  real<lower=0>          pri_p_die_if_sev_a;
  real<lower=0>          pri_p_die_if_sev_b;
  real<lower=0>          pri_p_die_if_sym_a;
  real<lower=0>          pri_p_die_if_sym_b;
  // p(diag)  
  real<lower=0>          pri_p_diag_if_sym_a;
  real<lower=0>          pri_p_diag_if_sym_b;
  real<lower=0>          pri_p_diag_if_sev_a;
  real<lower=0>          pri_p_diag_if_sev_b;
  // delay to progression
  real<lower=0>          inf_prg_delay_shap; //~~ new names
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
  real<lower=0>          pri_cas_rep_delay_shp_a;
  real<lower=0>          pri_cas_rep_delay_shp_b;
  real<lower=0>          pri_die_rep_delay_shp_a;
  real<lower=0>          pri_die_rep_delay_shp_b;


  // turn on/off negative binomial
  int<lower = 0, upper = 1> nb_yes; 
  int<lower = 0, upper = 1> obs_cas_rep;
  int<lower = 0, upper = 1> obs_die_rep;
}
///////////////////////////////////////////////////////////
transformed data {

  int<lower=1>  N_days_tot;
  
  N_days_tot = N_days + N_days_before; 

}
///////////////////////////////////////////////////////////
parameters {
  
// INCIDENCE  
 real                  log_new_inf_0; // intercept in log space
 vector[N_days_tot-1]  deriv1_log_new_inf; // first derivative of the random walk 

// SYMPTOMS AND CARE 
  real<lower=0, upper=1>  p_sym_if_inf;
  real<lower=0, upper=1>  p_sev_if_sym;
  real<lower=0, upper=1>  p_die_if_sev;
  real<lower=0, upper=1>  p_die_if_sym;

  real<lower=0, upper =1>     scale_dx_delay_sym; 
  real<lower=0, upper =1>     scale_dx_delay_sev; 
  
  real<lower=0>     cas_rep_delay_shap;
  real<lower=0>     cas_rep_delay_rate;
  real<lower=0>     die_rep_delay_shap;
  real<lower=0>     die_rep_delay_rate;

// DIAGNOSIS // probability of diagnosis at each illness state
  real<lower=0, upper=1>  p_diag_if_sym;
  real<lower=0, upper=1>  p_diag_if_sev;

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

  vector[N_days_tot-2]    deriv2_log_new_inf;
  
// delay probabilitily of reporting w x days delay (PDF of delay distribution)

  vector[N_days_tot]  inf_prg_delay;
  vector[N_days_tot]  sym_prg_delay;
  vector[N_days_tot]  sev_prg_delay;
  
  vector[N_days_tot]  sym_diag_delay;
  vector[N_days_tot]  sev_diag_delay;
  
  vector[N_days_tot]  cas_rep_delay;
  vector[N_days_tot]  die_rep_delay;
  
  // cumulative report delays
  vector[N_days_tot]  cas_cum_report_delay; 
  vector[N_days_tot]  die_cum_report_delay; 
  
  vector[N_days_tot]  inf_cum_prg_delay;
  vector[N_days_tot]  sym_cum_prg_delay; 
  vector[N_days_tot]  sev_cum_prg_delay; 
  
// OUTCOMES
  vector[N_days_tot]  new_sym; // new in state per day
  vector[N_days_tot]  new_sev;
  vector[N_days_tot]  new_die;
  
  vector[N_days_tot]  cumul_sym;
  vector[N_days_tot]  cumul_die; 

  vector[N_days_tot]  new_sym_dx; // enter new state << undiagnosed >>
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

// DLEAYS /////////////////////////

//  progression
  for(i in 1:N_days_tot) {
    inf_prg_delay[i] = gamma_cdf(i+0.0, inf_prg_delay_shap, inf_prg_delay_rate)
      - gamma_cdf(i-1.0, inf_prg_delay_shap, inf_prg_delay_rate);
    sym_prg_delay[i] = gamma_cdf(i+0.0, sym_prg_delay_shap, sym_prg_delay_rate)
      - gamma_cdf(i-1.0, sym_prg_delay_shap, sym_prg_delay_rate);
    sev_prg_delay[i] = gamma_cdf(i+0.0, sev_prg_delay_shap, sev_prg_delay_rate)
      - gamma_cdf(i-1.0, sev_prg_delay_shap, sev_prg_delay_rate);
  }
  
for(i in 1:N_days_tot){
    sym_diag_delay[i] = gamma_cdf(i+0.0, sym_prg_delay_shap, sym_prg_delay_rate/scale_dx_delay_sym)
      - gamma_cdf(i-1.0, sym_prg_delay_shap, sym_prg_delay_rate/scale_dx_delay_sym);
    sev_diag_delay[i] = gamma_cdf(i+0.0, sev_prg_delay_shap, sev_prg_delay_rate/scale_dx_delay_sev)
      - gamma_cdf(i-1.0, sev_prg_delay_shap, sev_prg_delay_rate/scale_dx_delay_sev);
}
  
// reporting delays//~~
  for(i in 1:N_days_tot) {
    cas_rep_delay[i] = gamma_cdf(i+0.0, cas_rep_delay_shap,cas_rep_delay_rate)
      - gamma_cdf(i-1.0, cas_rep_delay_shap, cas_rep_delay_rate);
    die_rep_delay[i] = gamma_cdf(i+0.0, die_rep_delay_shap, die_rep_delay_rate)
      - gamma_cdf(i-1.0, die_rep_delay_shap, die_rep_delay_rate);
}

// cumlative delays
cas_cum_report_delay = cumulative_sum(cas_rep_delay);
die_cum_report_delay = cumulative_sum(die_rep_delay);

inf_cum_prg_delay = cumulative_sum(inf_prg_delay);
sym_cum_prg_delay = cumulative_sum(sym_prg_delay);
sev_cum_prg_delay = cumulative_sum(sev_prg_delay);

// CASCADE OF INCIDENT OUTCOMES (TOTAL) ///////////////////////////////////
  new_sym = rep_vector(0, N_days_tot);
  new_sev = rep_vector(0, N_days_tot);
  new_die = rep_vector(0, N_days_tot);
  
  // INCIDENCE

  log_new_inf[1] = log_new_inf_0;
   for(i in 1:(N_days_tot-1)) {
    log_new_inf[i+1] =  log_new_inf[i] + deriv1_log_new_inf[i];
  }
  for(i in 1:(N_days_tot-2)) {
    deriv2_log_new_inf[i] = log_new_inf[i+1] * 2 - log_new_inf[i] - log_new_inf[i+2];
  }

  new_inf = exp(log_new_inf);
  
  for(i in 1:N_days_tot) {
    for(j in 1:(N_days_tot-i+1)) {
        new_sym[i+(j-1)] += new_inf[i] * p_sym_if_inf * inf_prg_delay[j];
      }
    }

  for(i in 1:N_days_tot) {
    for(j in 1:(N_days_tot-i+1)) {
        new_sev[i+(j-1)] += new_sym[i] * p_sev_if_sym * sym_prg_delay[j];
      }
    }

  for(i in 1:N_days_tot) {
    for(j in 1:(N_days_tot-i+1)) {
        new_die[i+(j-1)] += new_sev[i] * p_die_if_sev * sev_prg_delay[j];
      }
    }
cumul_sym = cumulative_sum(new_sym); 
cumul_die = cumulative_sum(new_die); 

cumul_die[N_days_tot] = cumul_sym[N_days_tot] * p_die_if_sym; 

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
    for(j in 1:(N_days_tot-i+1)){
      new_sym_dx[i+(j-1)] += new_sym[i]  * 
                             p_sev_if_sym * (1 - sym_cum_prg_delay[j]) *
                             p_diag_if_sym * 
                             sym_diag_delay[j];
    }
  }
  
    // cascade from diagnosis 
          for(i in 1:N_days_tot){
            for(j in 1:(N_days_tot-i+1)){
              dx_sym_sev[i+(j-1)] += new_sym_dx[i] * p_sev_if_sym *
              sym_prg_delay[j];
            }
          }
                  for(i in 1:N_days_tot){
            for(j in 1:(N_days_tot-i+1)){
              dx_sym_die[i+(j-1)] += dx_sym_sev[i] * p_die_if_sev *
                sev_prg_delay[j];
            }
          }    
    
  // diagnosed at severe 
  for(i in 1:N_days_tot){
    for(j in 1:(N_days_tot-i+1)){
      new_sev_dx[i+(j-1)] += (new_sev[i] - dx_sym_sev[i]) * 
                             p_die_if_sev * (1 - sev_cum_prg_delay[j]) * 
                             p_diag_if_sev * 
                             sev_diag_delay[j]; 
    }
  }
  
  for(i in 1:N_days_tot){
    for(j in 1:(N_days_tot-i+1)){
      dx_sev_die[i+(j-1)] += new_sev_dx[i] * p_die_if_sev * sev_prg_delay[j];
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

for(i in 1:N_days_tot)  {
  occur_cas[i] += diag_all[i] * cas_cum_report_delay[N_days_tot - i + 1];
  occur_die[i] += new_die_dx[i] * die_cum_report_delay[N_days_tot - i + 1];
}

}
///////////////////////////////////////////////////////////  
model {
//// PRIORS
  // INCIDENCE
  log_new_inf_0         ~ normal(pri_log_new_inf_0_mu, pri_log_new_inf_0_sd); // 
  deriv1_log_new_inf    ~ normal(0, pri_deriv1_log_new_inf_sd); //~~
  deriv2_log_new_inf    ~ normal(0, pri_deriv2_log_new_inf_sd);
  
  // SYMPTOMS AND CARE
  p_sym_if_inf         ~ beta(pri_p_sym_if_inf_a, pri_p_sym_if_inf_b);
  p_sev_if_sym         ~ beta(pri_p_sev_if_sym_a, pri_p_sev_if_sym_b);
  p_die_if_sev         ~ beta(pri_p_die_if_sev_a, pri_p_die_if_sev_b);
  p_die_if_sym         ~ beta(pri_p_die_if_sym_a, pri_p_die_if_sym_b);

  scale_dx_delay_sym   ~ beta(scale_dx_delay_sym_a, scale_dx_delay_sym_b); 
  scale_dx_delay_sev   ~ beta(scale_dx_delay_sev_a, scale_dx_delay_sev_b);
  
  cas_rep_delay_shap   ~ lognormal(log(pri_cas_rep_delay_shp_a), 0.5);  
  cas_rep_delay_rate   ~ lognormal(log(pri_cas_rep_delay_shp_b), 0.5); 
  die_rep_delay_shap   ~ lognormal(log(pri_die_rep_delay_shp_a), 0.5);  
  die_rep_delay_rate   ~ lognormal(log(pri_die_rep_delay_shp_b), 0.5);

// DIAGNOSIS    
  p_diag_if_sym        ~ beta(pri_p_diag_if_sym_a, pri_p_diag_if_sym_b);
  p_diag_if_sev        ~ beta(pri_p_diag_if_sev_a, pri_p_diag_if_sev_b);
  
  inv_sqrt_phi_c       ~ normal(0, 1); //~~
  inv_sqrt_phi_d       ~ normal(0, 1); //~~
    
//// LIKELIHOOD
//SWITCH TO NEG BIN
  if (nb_yes == 1) {
         for(i in 1:N_days) {
            obs_cas[i] ~ neg_binomial_2(occur_cas[i + N_days_before], phi_cas);
          }
       for(i in 1:N_days) {
          obs_die[i] ~ neg_binomial_2(occur_die[i + N_days_before], phi_die);
         }
// eventually you'll want cumulative cases
  } else { 
       for(i in 1:N_days) {
            obs_cas[i] ~ poisson(occur_cas[i + N_days_before]);
         }
 
       for(i in 1:N_days) {
          obs_die[i] ~ poisson(occur_die[i + N_days_before]);
          }

   }
}
///////////////////////////////////////////////////////////
generated quantities {
}
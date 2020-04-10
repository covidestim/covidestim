
################# STILL TO DO #################

# -- change in testing probability by time
# -- age
# -- syndromic surveillance inputs
# -- incorporate data on test negatives --> indicator of testing access?


stan_code <- "

///////////////////////////////////////////////////////////
data {
///~~~~~~~ Define ~~~~~~~
int<lower=0>     N_conf_cases;
int<lower=0>     N_days;
int<lower=0>     N_days_extra;
int<lower=0>     Max_delay;
int<lower=0>     cases_test_day[N_conf_cases]; // test date
int<lower=0>     cases_days_delay[N_conf_cases]; // delay to diag

real              pri_log_new_inf_0_mu;
real<lower=0>     pri_log_new_inf_0_sd;
real              pri_log_new_inf_drift_mu;
real<lower=0>     pri_log_new_inf_drift_sd;
real<lower=0>     pri_sigma_deriv1_log_new_inf_sd;
real<lower=0>     pri_deriv2_log_new_inf_sd;
// priors on conditional model state
real<lower=0>     pri_p_sym_if_inf_mn;
real<lower=0>     pri_p_sym_if_inf_ss;
real<lower=0>     pri_p_hos_if_sym_mn;
real<lower=0>     pri_p_hos_if_sym_ss;
real<lower=0>     pri_p_die_if_hos_mn;
real<lower=0>     pri_p_die_if_hos_ss;

// priors on delays in progression to next model state
real<lower=0>     pri_inf_prg_delay_mn_mn; // simplify
real<lower=0>     pri_inf_prg_delay_mn_cv; // simplify
real<lower=0>     inf_prg_delay_cv;

real<lower=0>     pri_sym_prg_delay_mn_mn; // simplify
real<lower=0>     pri_sym_prg_delay_mn_cv; // simplify
real<lower=0>     sym_prg_delay_cv;

real<lower=0>     pri_hos_prg_delay_mn_mn; // simplify
real<lower=0>     pri_hos_prg_delay_mn_cv; // simplify
real<lower=0>     hos_prg_delay_cv;

// priors on resoultion delay (state -> recovered)
real<lower=0>     pri_inf_res_delay_mn_mn; // simplify
real<lower=0>     pri_inf_res_delay_mn_cv; // simplify
real<lower=0>     inf_res_delay_cv;

real<lower=0>     pri_sym_res_delay_mn_mn; // simplify
real<lower=0>     pri_sym_res_delay_mn_cv; // simplify
real<lower=0>     sym_res_delay_cv;

real<lower=0>     pri_hos_res_delay_mn_mn; // simplify
real<lower=0>     pri_hos_res_delay_mn_cv; // simplify
real<lower=0>     hos_res_delay_cv;

// reporting delays
real<lower=0>     pri_report_delay_mn_mn; // simplify
real<lower=0>     pri_report_delay_mn_cv; // simplify

real<lower=0>     pri_report_delay_cv_mn; // simplify
real<lower=0>     pri_report_delay_cv_cv; // simplify

// diagnosis delays
real<lower=0>     pri_p_diag_if_inf_mn;
real<lower=0>     pri_p_diag_if_inf_ss;

real<lower=0>     pri_p_diag_if_sym_mn;
real<lower=0>     pri_p_diag_if_sym_ss;

real<lower=0>     pri_p_diag_if_hos_mn;
real<lower=0>     pri_p_diag_if_hos_ss;

int<lower = 0, upper = 1> nb_yes; // turns on a negative binomial (currently fit poisson)

}
///////////////////////////////////////////////////////////
transformed data {
///~~~~~~~ Define ~~~~~~~
  int<lower=1>  N_days_tot;
  int           rep_tri_conf_cases[N_days, Max_delay+1]; // reporting triangle -- laying out cases by date of report and by delay. Rows are date of report, columns are delay.
  ///////////////////// PRIORS
  // p(progression)
  real          pri_p_sym_if_inf_a; // all these (through 136) take distributions entered as data, turns them into alpha, beta parameters of gamma distribution
  real          pri_p_sym_if_inf_b;
  
  real          pri_p_hos_if_sym_a;
  real          pri_p_hos_if_sym_b;
  
  real          pri_p_die_if_hos_a;
  real          pri_p_die_if_hos_b;
// time to progression
  real          pri_inf_prg_delay_mn_shap;
  real          pri_inf_prg_delay_mn_rate;
  
  real          pri_sym_prg_delay_mn_shap;
  real          pri_sym_prg_delay_mn_rate;
  
  real          pri_hos_prg_delay_mn_shap;
  real          pri_hos_prg_delay_mn_rate;
// time to recovered 
  real          pri_inf_res_delay_mn_shap;
  real          pri_inf_res_delay_mn_rate;
  
  real          pri_sym_res_delay_mn_shap;
  real          pri_sym_res_delay_mn_rate;
  
  real          pri_hos_res_delay_mn_shap;
  real          pri_hos_res_delay_mn_rate;
// report delay
// shape and mean of delay distribution
// gamma priors on mean and cv of delay distribution
  real          pri_report_delay_mn_shap;
  real          pri_report_delay_mn_rate;
  real          pri_report_delay_cv_shap;
  real          pri_report_delay_cv_rate;
// p(diag)  
  real          pri_p_diag_if_inf_a;
  real          pri_p_diag_if_inf_b;
  
  real          pri_p_diag_if_sym_a;
  real          pri_p_diag_if_sym_b;
  
  real          pri_p_diag_if_hos_a;
  real          pri_p_diag_if_hos_b;
// delay to progression
  real          inf_prg_delay_shap;
  real          sym_prg_delay_shap;
  real          hos_prg_delay_shap;
// delay to recovered  
  real          inf_res_delay_shap;
  real          sym_res_delay_shap;
  real          hos_res_delay_shap;

///~~~~~~~ Assign values ~~~~~~~
  pri_p_sym_if_inf_a = pri_p_sym_if_inf_mn*pri_p_sym_if_inf_ss;
  pri_p_sym_if_inf_b = (1-pri_p_sym_if_inf_mn)*pri_p_sym_if_inf_ss;
  
  pri_p_hos_if_sym_a = pri_p_hos_if_sym_mn*pri_p_hos_if_sym_ss;
  pri_p_hos_if_sym_b = (1-pri_p_hos_if_sym_mn)*pri_p_hos_if_sym_ss;
  
  pri_p_die_if_hos_a = pri_p_die_if_hos_mn*pri_p_die_if_hos_ss;
  pri_p_die_if_hos_b = (1-pri_p_die_if_hos_mn)*pri_p_die_if_hos_ss;

  pri_inf_prg_delay_mn_shap = pow(pri_inf_prg_delay_mn_cv,-2); 
  pri_inf_prg_delay_mn_rate = pow(pri_inf_prg_delay_mn_cv,-2)/pri_inf_prg_delay_mn_mn;
  
  pri_sym_prg_delay_mn_shap = pow(pri_sym_prg_delay_mn_cv,-2);
  pri_sym_prg_delay_mn_rate = pow(pri_sym_prg_delay_mn_cv,-2)/pri_sym_prg_delay_mn_mn;
  
  pri_hos_prg_delay_mn_shap = pow(pri_hos_prg_delay_mn_cv,-2);
  pri_hos_prg_delay_mn_rate = pow(pri_hos_prg_delay_mn_cv,-2)/pri_hos_prg_delay_mn_mn;
 
  pri_inf_res_delay_mn_shap = pow(pri_inf_res_delay_mn_cv,-2);
  pri_inf_res_delay_mn_rate = pow(pri_inf_res_delay_mn_cv,-2)/pri_inf_res_delay_mn_mn;
  
  pri_sym_res_delay_mn_shap = pow(pri_sym_res_delay_mn_cv,-2);
  pri_sym_res_delay_mn_rate = pow(pri_sym_res_delay_mn_cv,-2)/pri_sym_res_delay_mn_mn;
  
  pri_hos_res_delay_mn_shap = pow(pri_hos_res_delay_mn_cv,-2);
  pri_hos_res_delay_mn_rate = pow(pri_hos_res_delay_mn_cv,-2)/pri_hos_res_delay_mn_mn;

  pri_report_delay_mn_shap = pow(pri_report_delay_mn_cv,-2);
  pri_report_delay_mn_rate = pow(pri_report_delay_mn_cv,-2)/pri_report_delay_mn_mn;
  pri_report_delay_cv_shap = pow(pri_report_delay_cv_cv,-2);
  pri_report_delay_cv_rate = pow(pri_report_delay_cv_cv,-2)/pri_report_delay_cv_mn;
  
  pri_p_diag_if_inf_a = pri_p_diag_if_inf_mn*pri_p_diag_if_inf_ss;
  pri_p_diag_if_inf_b = (1-pri_p_diag_if_inf_mn)*pri_p_diag_if_inf_ss;
  
  pri_p_diag_if_sym_a = pri_p_diag_if_sym_mn*pri_p_diag_if_sym_ss;
  pri_p_diag_if_sym_b = (1-pri_p_diag_if_sym_mn)*pri_p_diag_if_sym_ss;
  
  pri_p_diag_if_hos_a = pri_p_diag_if_hos_mn*pri_p_diag_if_hos_ss;
  pri_p_diag_if_hos_b = (1-pri_p_diag_if_hos_mn)*pri_p_diag_if_hos_ss;

  inf_prg_delay_shap = pow(inf_prg_delay_cv,-2);
  sym_prg_delay_shap = pow(sym_prg_delay_cv,-2);
  hos_prg_delay_shap = pow(hos_prg_delay_cv,-2);

  inf_res_delay_shap = pow(inf_res_delay_cv,-2);
  sym_res_delay_shap = pow(sym_res_delay_cv,-2);
  hos_res_delay_shap = pow(hos_res_delay_cv,-2);
  
//////////////////// End Priors   

  N_days_tot = N_days + N_days_extra; // how many days to project for

  for(i in 1:N_days){
    for(j in 1:(Max_delay+1)){
      rep_tri_conf_cases[i,j] = 0;
    }
  }

  for(i in 1:N_conf_cases){
    rep_tri_conf_cases[cases_test_day[i], cases_days_delay[i]+1] += 1;
  }

}
///////////////////////////////////////////////////////////
parameters {
///~~~~~~~ Define ~~~~~~~
// INCIDENCE // building the random walk 
  real                  log_new_inf_0; // intercept in log space
  real                  log_new_inf_drift; // mean day on day change
  vector[N_days_tot-1]  deriv1_log_new_inf; // first derivative of the random walk 
  real<lower=0>         sigma_deriv1_log_new_inf; // parameter for the SD of the rw
//  real<lower=0> sigma_deriv2_log_new_inf;

// SYMPTOMS AND CARE // probability of transitioning to next illnesses state
// currently assume a fixed cv
  real<lower=0, upper=1> p_sym_if_inf;
  real<lower=0, upper=1> p_hos_if_sym;
  real<lower=0, upper=1> p_die_if_hos;
// delay associated with transition to next illness state
  real<lower=0> inf_prg_delay_mn;
  real<lower=0> sym_prg_delay_mn;
  real<lower=0> hos_prg_delay_mn;
// delay associated with transition to recovered 
  real<lower=0> inf_res_delay_mn;
  real<lower=0> sym_res_delay_mn;
  real<lower=0> hos_res_delay_mn;

// DIAGNOSIS // probability of diagnosis at each illness state
  real<lower=0, upper=1> p_diag_if_inf;
  real<lower=0, upper=1> p_diag_if_sym;
  real<lower=0, upper=1> p_diag_if_hos;
// REPORTING DElAYS  
  real<lower=0> report_delay_mn;
  real<lower=0> report_delay_cv;

// LIKELIHOOD
  real<lower=0> inv_sqrt_phi;

}
///////////////////////////////////////////////////////////
transformed parameters {
///~~~~~~~ Define ~~~~~~~
// INCIDENCE
  vector[N_days_tot]     log_new_inf;
  vector[N_days_tot]     new_inf;
  vector[N_days_tot-2]   deriv2_log_new_inf;

// DELAYS probabilitily of reporting w x days delay (PDF of delay distribution)
  vector[N_days_tot]  inf_prg_delay;
  vector[N_days_tot]  sym_prg_delay;
  vector[N_days_tot]  hos_prg_delay;

  vector[N_days_tot]  inf_res_delay;
  vector[N_days_tot]  sym_res_delay;
  vector[N_days_tot]  hos_res_delay;

  vector[N_days_tot]  report_delay;
// rate parameter -- can only be estimated from mean (the parameter)
// parameter of delay distribution
  real<lower=0> inf_prg_delay_rate;
  real<lower=0> sym_prg_delay_rate;
  real<lower=0> hos_prg_delay_rate;

  real<lower=0> inf_res_delay_rate;
  real<lower=0> sym_res_delay_rate;
  real<lower=0> hos_res_delay_rate;

  real<lower=0> report_delay_shap;
  real<lower=0> report_delay_rate;

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

  vector[N_days_tot]  cur_inf_u; // currently in state << undiagnosed >>
  vector[N_days_tot]  cur_sym_u;
  vector[N_days_tot]  cur_hos_u;

  vector[N_days_tot]  cum_diag_inf; // cumulative diagnoses by state
  vector[N_days_tot]  cum_diag_sym;
  vector[N_days_tot]  cum_diag_hos;

  vector[N_days_tot]  diag_inf; // incident diagnosis
  vector[N_days_tot]  diag_sym;
  vector[N_days_tot]  diag_hos;
  vector[N_days_tot]  diag_all;

  real                phi;

  matrix[N_days,Max_delay+1]  rep_tri_conf_cases_mu;

///~~~~~~~ Assign values ~~~~~~~

  phi = pow(inv_sqrt_phi,-2);

// INCIDENCE
  log_new_inf[1] = log_new_inf_0;
  for(i in 1:(N_days_tot-1)) {
  //  log_new_inf[i+1] =  log_new_inf[i] + log_new_inf_drift + deriv1_log_new_inf[i] * sigma_deriv1_log_new_inf;
    log_new_inf[i+1] =  log_new_inf[i] + log_new_inf_drift + deriv1_log_new_inf[i] ;
  }

  new_inf = exp(log_new_inf);

  for(i in 1:(N_days_tot-2)) {
    deriv2_log_new_inf[i] = log_new_inf[i+1] * 2 - log_new_inf[i] - log_new_inf[i+2];
  }

// DELAYS // progression
  inf_prg_delay_rate = inf_prg_delay_shap/inf_prg_delay_mn;
  sym_prg_delay_rate = sym_prg_delay_shap/sym_prg_delay_mn;
  hos_prg_delay_rate = hos_prg_delay_shap/hos_prg_delay_mn;

  for(i in 1:N_days_tot) {
    inf_prg_delay[i] = gamma_cdf(i+0.0, inf_prg_delay_shap, inf_prg_delay_rate) - gamma_cdf(i-1.0, inf_prg_delay_shap, inf_prg_delay_rate);
    sym_prg_delay[i] = gamma_cdf(i+0.0, sym_prg_delay_shap, sym_prg_delay_rate) - gamma_cdf(i-1.0, sym_prg_delay_shap, sym_prg_delay_rate);
    hos_prg_delay[i] = gamma_cdf(i+0.0, hos_prg_delay_shap, hos_prg_delay_rate) - gamma_cdf(i-1.0, hos_prg_delay_shap, hos_prg_delay_rate);
  }
// DELAYS // resolution of case
  inf_res_delay_rate = inf_res_delay_shap/inf_res_delay_mn;
  sym_res_delay_rate = sym_res_delay_shap/sym_res_delay_mn;
  hos_res_delay_rate = hos_res_delay_shap/hos_res_delay_mn;

  for(i in 1:N_days_tot) {
    inf_res_delay[i] = gamma_cdf(i+0.0, inf_res_delay_shap, inf_res_delay_rate) - gamma_cdf(i-1.0, inf_res_delay_shap, inf_res_delay_rate);
    sym_res_delay[i] = gamma_cdf(i+0.0, sym_res_delay_shap, sym_res_delay_rate) - gamma_cdf(i-1.0, sym_res_delay_shap, sym_res_delay_rate);
    hos_res_delay[i] = gamma_cdf(i+0.0, hos_res_delay_shap, hos_res_delay_rate) - gamma_cdf(i-1.0, hos_res_delay_shap, hos_res_delay_rate);
  }
// DELAYS // reporting
  report_delay_shap  = pow(report_delay_cv,-2);
  report_delay_rate  = pow(report_delay_cv,-2)/report_delay_mn;

  for(i in 1:N_days_tot) {
    report_delay[i] = gamma_cdf(i+0.0, report_delay_shap, report_delay_rate) - gamma_cdf(i-1.0, report_delay_shap, report_delay_rate);
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
  // could be "prof_inf" <- infections that progress
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
        res_inf_u[i+(j-1)] += new_inf[i] * (1-p_sym_if_inf) * inf_res_delay[j] * pow(1-p_diag_if_inf,j-1);
      }
    }
  }

  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_sym_u[i+(j-1)] += new_inf[i] * p_sym_if_inf * inf_prg_delay[j] * pow(1-p_diag_if_inf,j-1);
      }
    }
  }

  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        res_sym_u[i+(j-1)] += new_sym_u[i] * (1-p_hos_if_sym) * sym_res_delay[j] * pow(1-p_diag_if_sym,j-1);
      }
    }
  }

  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_hos_u[i+(j-1)] += new_sym_u[i] * p_hos_if_sym * sym_prg_delay[j] * pow(1-p_diag_if_sym,j-1);
      }
    }
  }

  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        res_hos_u[i+(j-1)] += new_hos_u[i] * (1-p_die_if_hos) * hos_res_delay[j] * pow(1-p_diag_if_hos,j-1);
      }
    }
  }


  for(i in 1:N_days_tot) {
    for(j in 1:N_days_tot) {
      if(i+(j-1) <= N_days_tot){
        new_die_u[i+(j-1)] += new_hos_u[i] * p_die_if_hos * hos_prg_delay[j] * pow(1-p_diag_if_hos,j-1);
      }
    }
  }

// CURRENT STATUS  << UNDIAGNOSED >>
  cur_inf_u = rep_vector(0, N_days_tot);
  cur_sym_u = rep_vector(0, N_days_tot);
  cur_hos_u = rep_vector(0, N_days_tot);

  cur_inf_u[1] = new_inf[1] - new_sym_u[1] - res_inf_u[1];
  for(i in 2:N_days_tot) {
    cur_inf_u[i] = cur_inf_u[i-1]*(1-p_diag_if_inf) + new_inf[i] - new_sym_u[i] - res_inf_u[i];
  }

  cur_sym_u[1] = new_sym_u[1] - new_hos_u[1] - res_sym_u[1];
  for(i in 2:N_days_tot) {
    cur_sym_u[i] = cur_sym_u[i-1]*(1-p_diag_if_sym) + new_sym_u[i] - new_hos_u[i] - res_sym_u[i]; 
  }

  cur_hos_u[1] = new_hos_u[1] - new_die_u[1] - res_hos_u[1];  
  for(i in 2:N_days_tot) {
    cur_hos_u[i] = cur_hos_u[i-1]*(1-p_diag_if_hos) + new_hos_u[i] - new_die_u[i] - res_hos_u[i];
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
///////////////////////////////////////////////////////////  
model {
///~~~~~~~ Assign values ~~~~~~~
//// PRIORS
// INCIDENCE
  log_new_inf_0             ~ normal(pri_log_new_inf_0_mu, pri_log_new_inf_0_sd); // 
  log_new_inf_drift         ~ normal(0, pri_log_new_inf_drift_sd); // 
  sigma_deriv1_log_new_inf  ~ normal(0, pri_sigma_deriv1_log_new_inf_sd);  // should be cauchy // sd of randomw walk
  deriv1_log_new_inf        ~ normal(0, 1);
  deriv2_log_new_inf        ~ normal(0, pri_deriv2_log_new_inf_sd);

// SYMPTOMS AND CARE
  p_sym_if_inf              ~ beta(pri_p_sym_if_inf_a, pri_p_sym_if_inf_b);
  p_hos_if_sym              ~ beta(pri_p_hos_if_sym_a, pri_p_hos_if_sym_b);
  p_die_if_hos              ~ beta(pri_p_die_if_hos_a, pri_p_die_if_hos_b);

  inf_prg_delay_mn          ~ gamma(pri_inf_prg_delay_mn_shap, pri_inf_prg_delay_mn_rate);
  sym_prg_delay_mn          ~ gamma(pri_sym_prg_delay_mn_shap, pri_sym_prg_delay_mn_rate);
  hos_prg_delay_mn          ~ gamma(pri_hos_prg_delay_mn_shap, pri_hos_prg_delay_mn_rate);

  inf_res_delay_mn          ~ gamma(pri_inf_res_delay_mn_shap, pri_inf_res_delay_mn_rate);
  sym_res_delay_mn          ~ gamma(pri_sym_res_delay_mn_shap, pri_sym_res_delay_mn_rate);
  hos_res_delay_mn          ~ gamma(pri_hos_res_delay_mn_shap, pri_hos_res_delay_mn_rate);

  report_delay_mn           ~ gamma(pri_report_delay_mn_shap, pri_report_delay_mn_rate); 
  report_delay_cv           ~ gamma(pri_report_delay_cv_shap, pri_report_delay_cv_rate);

  p_diag_if_inf              ~ beta(pri_p_diag_if_inf_a, pri_p_diag_if_inf_b);
  p_diag_if_sym              ~ beta(pri_p_diag_if_sym_a, pri_p_diag_if_sym_b);
  p_diag_if_hos              ~ beta(pri_p_diag_if_hos_a, pri_p_diag_if_hos_b);

  inv_sqrt_phi               ~ normal(0, 1); 

////   LIKELIHOOD
// REPORTED CASES
  for(i in 1:N_days) {
    for(j in 1:(Max_delay+1)) {
      if((i+j) < (N_days+2)){
        if(nb_yes==1){
          rep_tri_conf_cases[i,j] ~ neg_binomial_2(rep_tri_conf_cases_mu[i,j], phi);
        } else {
          rep_tri_conf_cases[i,j] ~ poisson(rep_tri_conf_cases_mu[i,j]);
        }
      }
    }
  }

}
///////////////////////////////////////////////////////////
generated quantities {

}

"

###############################
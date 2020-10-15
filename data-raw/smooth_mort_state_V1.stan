stan_code <- "
data {
  // INPUT DATA
  int<lower=0>           N_obs; 
  int<lower=0>           N_state; 
  int<lower=0>           N_age; 
  int<lower=0>           obs_die[N_obs]; // vector of deaths
  int<lower=0>           state[N_obs]; // vector of state factors
  int<lower=0>           age[N_obs]; // vector of age factors
  int<lower=0>           age_state[N_obs]; // vector of state*age factors
  real<lower=0>          sigma_re_sd;   

}
///////////////////////////////////////////////////////////
transformed data {
}
///////////////////////////////////////////////////////////
parameters {
// fixed effects 
  real            b_0; //  intercept
  real            b_age0[N_age-1]; // fe for age main effect
  real            b_state0[N_state-1]; // fe for state main effect
  // random effects 
  real<lower=0>   sig_state_age_re; //  sigma for re
  real            b_state_age_re[N_age*N_state]; // re for state*age
}
///////////////////////////////////////////
transformed parameters {
///~~~~~~~ Define ~~~~~~~
  real            b_age[N_age]; // fe for age main effect
  real            b_state[N_state]; // fe for state main effect  vector[N_days_tot]      log_new_inf;
  real            pred[N_obs];
  
///~~~~~~~ Calc ~~~~~~~
  b_age[1] = 0;
  b_age[2:N_age] = b_age0;
  b_state[1] = 0;
  b_state[2:N_state] = b_state0;
  for(i in 1:N_obs){
    pred[i] = exp(b_0 + b_age[age[i]] + b_state[state[i]] + b_state_age_re[age_state[i]]*sig_state_age_re);
  }

}
///////////////////////////////////////////////////////////  
model {
///~~~~~~~ prior ~~~~~~~
// FE
  b_0              ~ normal(0, 10);
  b_age0           ~ normal(0, 10);
  b_state0         ~ normal(0, 10);
// RE
  b_state_age_re   ~ normal(0, 1);
  sig_state_age_re ~ normal(0, sigma_re_sd);
  
///~~~~~~~ likelihood ~~~~~~~
 obs_die ~ poisson(pred);
 
}
///////////////////////////////////////////////////////////
generated quantities {
}
"


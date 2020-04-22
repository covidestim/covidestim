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


///////////////////////////////////////////////////////////  
model {
//// PRIORS

// INCIDENCE
//if(rw_yes == 1){
  log_new_inf_0             ~ normal(pri_log_new_inf_0_mu, pri_log_new_inf_0_sd); // 
//  log_new_inf_drift         ~ normal(0, pri_log_new_inf_drift_sd); // 

  // should be cauchy // sd of randomw walk
  sigma_deriv1_log_new_inf  ~ normal(0, pri_sigma_deriv1_log_new_inf_sd);  
  deriv1_log_new_inf        ~ normal(0, 1);
  deriv2_log_new_inf        ~ normal(0, pri_deriv2_log_new_inf_sd);
//} else {
// deriv2_b_spline            ~ student_t(10, 0, 5);
// deriv1_b_spline            ~ student_t(10, 0, 1);
//}

// SYMPTOMS AND CARE
  p_sym_if_inf              ~ beta(pri_p_sym_if_inf_a, pri_p_sym_if_inf_b);
  p_hos_if_sym              ~ beta(pri_p_hos_if_sym_a, pri_p_hos_if_sym_b);
  p_die_if_hos              ~ beta(pri_p_die_if_hos_a, pri_p_die_if_hos_b);
  inf_prg_delay_mn          ~ gamma(inf_prg_delay_shap, inf_prg_delay_rate);
  sym_prg_delay_mn          ~ gamma(sym_prg_delay_shap, sym_prg_delay_rate);
  hos_prg_delay_mn          ~ gamma(hos_prg_delay_shap, hos_prg_delay_rate);
  inf_res_delay_mn          ~ gamma(inf_res_delay_shap, inf_res_delay_rate);
  sym_res_delay_mn          ~ gamma(sym_res_delay_shap, sym_res_delay_rate);
  hos_res_delay_mn          ~ gamma(hos_res_delay_shap, hos_res_delay_rate);
  report_delay_mn           ~ gamma(pri_report_delay_shap, pri_report_delay_rate); 
  p_diag_if_inf             ~ beta(pri_p_diag_if_inf_a, pri_p_diag_if_inf_b);
  p_diag_if_sym             ~ beta(pri_p_diag_if_sym_a, pri_p_diag_if_sym_b);
  p_diag_if_hos             ~ beta(pri_p_diag_if_hos_a, pri_p_diag_if_hos_b);
  //inv_sqrt_phi              ~ normal(0, 1); 
  
////   LIKELIHOOD
// REPORTED CASES
  for(i in 1:N_days) {
    for(j in 1:(Max_delay+1)) {
      if((i+j) < (N_days+2)){
       // if(nb_yes==1){
        //  rep_tri_conf_cases[i,j] ~ neg_binomial_2(rep_tri_conf_cases_mu[i,j], phi);
        //} else {
          rep_tri_conf_cases[i,j] ~ poisson(rep_tri_conf_cases_mu[i,j]);
        //}
      }
    }
  }
}

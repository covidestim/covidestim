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

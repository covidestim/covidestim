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

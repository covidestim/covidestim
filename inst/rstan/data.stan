  
///////////////////////////////////////////////////////////

data {
///~~~~~~~ Define ~~~~~~~
int<lower=0>     N_conf_cases;
int<lower=0>     N_days;
int<lower=0>     N_days_extra;
int<lower=0>     Max_delay;
int<lower=0>     cases_test_day[N_conf_cases]; // test date
int<lower=0>     cases_days_delay[N_conf_cases]; // delay to diag

real<lower=0>     pri_log_new_inf_0_mu;
real<lower=0>     pri_log_new_inf_0_sd;
real<lower=0>     pri_sigma_deriv1_log_new_inf_sd;
real<lower=0>     pri_deriv2_log_new_inf_sd;

//int<lower=1>     n_spl_par; // term for the spline
//matrix[(N_days+N_days_extra),n_spl_par] spl_basis; // term for the spline
//
//

// PRIORS
// p(progression)
real<lower=0>          pri_p_sym_if_inf_a; 
real<lower=0>          pri_p_sym_if_inf_b;
real<lower=0>          pri_p_hos_if_sym_a;
real<lower=0>          pri_p_hos_if_sym_b;
real<lower=0>          pri_p_die_if_hos_a;
real<lower=0>          pri_p_die_if_hos_b;
// delay to progression
real<lower=0>          inf_prg_delay_shap;
real<lower=0>          inf_prg_delay_rate;
real<lower=0>          sym_prg_delay_shap;
real<lower=0>          sym_prg_delay_rate;
real<lower=0>          hos_prg_delay_shap;
real<lower=0>          hos_prg_delay_rate;
// delay to recovered  
real<lower=0>          inf_res_delay_shap;
real<lower=0>          inf_res_delay_rate;
real<lower=0>          sym_res_delay_shap;
real<lower=0>          sym_res_delay_rate;
real<lower=0>          hos_res_delay_shap;
real<lower=0>          hos_res_delay_rate;
// report delay // to be simplified later
real<lower=0>          pri_report_delay_shap;
real<lower=0>          pri_report_delay_rate;
// p(diag)  
real<lower=0>          pri_p_diag_if_inf_a;
real<lower=0>          pri_p_diag_if_inf_b;
real<lower=0>          pri_p_diag_if_sym_a;
real<lower=0>          pri_p_diag_if_sym_b;
real<lower=0>          pri_p_diag_if_hos_a;
real<lower=0>          pri_p_diag_if_hos_b;

//int<lower = 0, upper = 1> nb_yes; // turns on a negative binomial (currently fit poisson)
//int<lower = 0, upper = 1> rw_yes; 
}

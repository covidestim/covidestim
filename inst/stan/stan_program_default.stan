functions {
  vector conv1d(vector x, vector kernel) {
    int nk = rows(kernel);
    int nx = rows(x);
    matrix[nx, nk] X;

    if (nx < nk)
      reject("nrow(x) must be >= nrow(kernel). x had nrow =", nx);


    // The first `nK-1` rows of the convolution are special, because we don't
    // yet have a full "window" of X entries to take the dot product of against
    // the kernel. We want to handle this case by taking the dot prod of
    // `X[1:i]` and the the last `i` elements in the kernel.
    //
    // By being clever about how we position the elements of `Xmatrix[1:i]`
    // within row `i` of `Y`, we can guarantee that during the final matmul
    // step (`Xmatrix * kernel`), the elements will be lined up with the
    // correct subsequence of the kernel.
    for (i in 1:nk) {

      // Fill the beginning of the column with the correct number of `0`
      // elements.
      //
      // What we are shooting for here is for the first column to have `nK-1`
      // `0`s in it, and for the last column to have no `0`s in it.
      if (i < nk)
        X[1:(nk - i), i] = rep_vector(0, nk - i);

      // Populate the rest of each column with the entries of X, until we run
      // out of slots in the matrix (at the bottom of the column). The only
      // column that should have a full set of `X` entries is the last column,
      // the `nK`'th column.
      X[(nk - i + 1):nx, i] = x[1:(nx - nk + i)];
    }

    return X * kernel;
  }

}

data {
  // INPUT DATA
  // int<lower=0>           N_days; // days of data
  // int<lower=0>           N_days_before; // days before data to init epi model
  int<lower=0>           N_weeks; // weeks of data
  int<lower=0>           N_weeks_before; // weeks before data to init epi model
  int<lower=0>           Max_delay; // maximum days delay 
  
  // int<lower=0>           obs_cas[N_days]; // vector of cases
  // int<lower=0>           obs_die[N_days]; // vector of deaths
  int<lower=0>           obs_cas[N_weeks]; // vector of cases
  int<lower=0>           obs_die[N_weeks]; // vector of deaths
  // vector<lower=0>[N_days + N_days_before]        obs_boost; // vector of booster data
  vector<lower=0>[N_weeks + N_weeks_before]        obs_boost; // vector of booster data
  real<lower=0>          pop_size; // population size
  
  int<lower=0>           N_ifr_adj; // length of ifr_adjustment
  vector<lower=0>[N_ifr_adj] ifr_adj; // ifr_adjustment
  // vector<lower=0>[N_days+N_days_before] ifr_vac_adj; // ifr_vaccine_adjustment
  vector<lower=0>[N_weeks+N_weeks_before] ifr_vac_adj; // ifr_vaccine_adjustment
  real<lower=0>          pri_ifr_decl_OR_a; 
  real<lower=0>          pri_ifr_decl_OR_b;
  real<lower=0>          ifr_adj_fixed;
  
  real<lower=0>          infect_dist_rate;
  real<lower=0>          infect_dist_shap;
  real<lower=0>          seropos_dist_rate;
  real<lower=0>          seropos_dist_shap;

  // terms for splines
  // spline parameters and bases
  int<lower=0>                              N_spl_par_rt;
  // matrix[N_days+N_days_before,N_spl_par_rt] spl_basis_rt;
  matrix[N_weeks+N_weeks_before,N_spl_par_rt] spl_basis_rt;
  int<lower=0>                              N_spl_par_dx;
  // matrix[N_days+N_days_before,N_spl_par_dx] spl_basis_dx;
  matrix[N_weeks+N_weeks_before,N_spl_par_dx] spl_basis_dx;

  // fixed delay distributions
  // fixed delay distribtions. time from inf -> sym, sym -> sev, sev -> die
  real<lower=0>          inf_prg_delay_shap; 
  real<lower=0>          inf_prg_delay_rate;
  real<lower=0>          asy_rec_delay_shap;
  real<lower=0>          asy_rec_delay_rate;
  real<lower=0>          sym_prg_delay_shap;
  real<lower=0>          sym_prg_delay_rate;
  real<lower=0>          sev_prg_delay_shap;
  real<lower=0>          sev_prg_delay_rate;

  // fixed delay from diagnosis to report
  real<lower=0>          cas_rep_delay_shap;
  real<lower=0>          cas_rep_delay_rate;
  real<lower=0>          die_rep_delay_shap;
  real<lower=0>          die_rep_delay_rate;
  
  //// control nobs
  // whether to assume zero cases and deaths during warm-up
  int<lower = 0, upper = 1> pre_period_zero; 

  //  what data are included --  cases, deaths: 
  int<lower = 0, upper = 1> cas_yes; 
  int<lower = 0, upper = 1> die_yes; 

  //  how are the data dated -- report, occurrence: 
  int<lower = 0, upper = 1> obs_cas_rep;
  int<lower = 0, upper = 1> obs_die_rep;

  //  how many days should be used for the moving average in the likelihood 
  //  function? 
  // int N_days_av; 
  // int<lower = 1, upper = 10> N_days_av; 

  // is there a last obeserved deaths data day?
  // int<lower=0> lastDeathDate;
  int<lower=0> lastDeathWeek;

  // is there a last obeserved case data day?
  // int<lower=0> lastCaseDate;
  int<lower=0> lastCaseWeek;

  /////////
  // TERMS FOR PRIOR DISTRIBTUIONS
  // for new infections
  real          pri_log_new_inf_0_mu;
  real<lower=0> pri_log_new_inf_0_sd;
  real          pri_logRt_mu;   
  real<lower=0> pri_logRt_sd;   
  real<lower=0> pri_serial_i_shap; 
  real<lower=0> pri_serial_i_rate; 
  real<lower=0> pri_deriv1_spl_par_sd;
  real<lower=0> pri_deriv2_spl_par_sd;
  
  // probabilities of progression inf -> sym -> sev -> die
  real<lower=0> pri_p_sym_if_inf_a; 
  real<lower=0> pri_p_sym_if_inf_b;
  real<lower=0> pri_p_sev_if_sym_a;
  real<lower=0> pri_p_sev_if_sym_b;
  real<lower=0> pri_p_die_if_sev_a;
  real<lower=0> pri_p_die_if_sev_b;

  // overall case fatality rate
  real<lower=0> pri_p_die_if_inf_a;
  real<lower=0> pri_p_die_if_inf_b;

  // probabilities of diagnosis 
     // rate ratio, pr(dx) asymptomatic to symptomatic
  real<lower=0> pri_rr_diag_asy_vs_sym_a; 
  real<lower=0> pri_rr_diag_asy_vs_sym_b;

      // rate ratio, pr(dx) symptomatic to severe
  real<lower=0> pri_rr_diag_sym_vs_sev_a; 
  real<lower=0> pri_rr_diag_sym_vs_sev_b;

     // probability of diagnosis at severe 
  real<lower=0> pri_p_diag_if_sev_a;
  real<lower=0> pri_p_diag_if_sev_b;

  // delay to diagnosis assumed to be some fraction of progression delay
  // Beta prior distribtuion for that fraction 
  real<lower=0> scale_dx_delay_sym_a; 
  real<lower=0> scale_dx_delay_sym_b; 
  real<lower=0> scale_dx_delay_sev_a; 
  real<lower=0> scale_dx_delay_sev_b;

  
  // input for the number of days to put the Rt prior on
  // int<lower=0>  N_days_pri_Rt;
  // real<lower=0> sd_pri_Rt;

}
///////////////////////////////////////////////////////////
transformed data {
  // int  N_days_tot;
  int  N_weeks_tot;

  // Moving sums
  // int<lower=0>           obs_cas_mvs[N_days]; // vector of cases
  // int<lower=0>           obs_die_mvs[N_days]; // vector of deaths
  // int<lower=0>           nda0 = N_days_av - 1;
  int<lower=0>           obs_cas_mvs[N_weeks]; // vector of cases
  int<lower=0>           obs_die_mvs[N_weeks]; // vector of deaths
  // real pop_sus; 
  // Progression delays
  vector[Max_delay]  inf_prg_delay_rv;
  vector[Max_delay]  asy_rec_delay_rv; 
  vector[Max_delay]  sym_prg_delay_rv;
  vector[Max_delay]  sev_prg_delay_rv;
 
  // Reporting delays
  vector[Max_delay]  cas_rep_delay_rv;
  vector[Max_delay]  die_rep_delay_rv;
 
  // Cumulative reporting delays
  // vector[N_days + N_days_before]  cas_cum_report_delay_rv; 
  // vector[N_days + N_days_before]  die_cum_report_delay_rv; 
  vector[N_weeks + N_weeks_before]  cas_cum_report_delay_rv; 
  vector[N_weeks + N_weeks_before]  die_cum_report_delay_rv; 
  
 //   int  idx1[N_days + N_days_before];
 // int  idx2[N_days + N_days_before];
 // vector[N_days + N_days_before] idx3;
   int  idx1[N_weeks + N_weeks_before];
 int  idx2[N_weeks + N_weeks_before];
 vector[N_weeks + N_weeks_before] idx3;
 // real log_new_inf0calc = log(pop_size * .5);
  // create 'N_days_tot', which is days of data plus days to model before first 
  // case or death 
  // N_days_tot = N_days + N_days_before; 
  N_weeks_tot = N_weeks + N_weeks_before; 
  // pop_sus = pop_size *.8;
  // Indexes for convolutions
// for(i in 1:N_days_tot) {
for(i in 1:N_weeks_tot) {
  if(i-Max_delay>0){
    idx1[i] = i-Max_delay+1;
    idx2[i] = 1;
  } else {
    idx1[i] = 1;
    idx2[i] = Max_delay-i+1;
  }
  // if(i < (N_days_tot - 1)){ ## 
  // idx3[i] = N_days_tot-1-i;
  if(i < (N_weeks_tot - 1)){ ## 
  idx3[i] = N_weeks_tot-1-i;
  } else {
    idx3[i] = 1;
  }
}

  // compute the moving sums
  // for(i in 1:N_days) {
  //     if(i < N_days_av) {
      //   obs_cas_mvs[i] = 0;
      //   obs_die_mvs[i] = 0;
      // } else {
      //   obs_cas_mvs[i] = sum(obs_cas[(i - nda0) : i]);
      //   obs_die_mvs[i] = sum(obs_die[(i - nda0) : i]);
      // }
// }
  for(i in 1:N_weeks) {
    obs_cas_mvs[i] = obs_cas[i];
    obs_die_mvs[i] = obs_die[i];
  }
 
  // calculate the daily probability of transitioning to a new disease state
  // for days 1 to 60 after entering that state
  for(i in 1:Max_delay) {
    inf_prg_delay_rv[1+Max_delay-i] =
      gamma_cdf(i   , inf_prg_delay_shap, inf_prg_delay_rate) -
      gamma_cdf(i-1 , inf_prg_delay_shap, inf_prg_delay_rate);

    asy_rec_delay_rv[1+Max_delay-i] =
      gamma_cdf(i   ,  asy_rec_delay_shap, asy_rec_delay_rate*2) -
      gamma_cdf(i-1 , asy_rec_delay_shap, asy_rec_delay_rate*2); 

    // diagnosis happens, on average, midway through the infectious period
    // therefore, we multiple the rate parameter by 2
    sym_prg_delay_rv[1+Max_delay-i] =
      gamma_cdf(i   , sym_prg_delay_shap, sym_prg_delay_rate) -
      gamma_cdf(i-1 , sym_prg_delay_shap, sym_prg_delay_rate);

    sev_prg_delay_rv[1+Max_delay-i] =
      gamma_cdf(i   , sev_prg_delay_shap, sev_prg_delay_rate) -
      gamma_cdf(i-1 , sev_prg_delay_shap, sev_prg_delay_rate);
  }
  
  // Calcluate the probability of reporting for each day after diagnosis
  // for 1 to 60 post diagnosis. 
  for(i in 1:Max_delay) {
    cas_rep_delay_rv[1+Max_delay-i] = 
      gamma_cdf(i   , cas_rep_delay_shap, cas_rep_delay_rate) -
      gamma_cdf(i-1 , cas_rep_delay_shap, cas_rep_delay_rate);

    die_rep_delay_rv[1+Max_delay-i] =
      gamma_cdf(i   , die_rep_delay_shap, die_rep_delay_rate) -
      gamma_cdf(i-1 , die_rep_delay_shap, die_rep_delay_rate);
  }
  
   // Make sure sum to 1
  inf_prg_delay_rv = inf_prg_delay_rv/sum(inf_prg_delay_rv);
  asy_rec_delay_rv = asy_rec_delay_rv/sum(asy_rec_delay_rv);
  sym_prg_delay_rv = sym_prg_delay_rv/sum(sym_prg_delay_rv);
  sev_prg_delay_rv = sev_prg_delay_rv/sum(sev_prg_delay_rv);
  cas_rep_delay_rv = cas_rep_delay_rv/sum(cas_rep_delay_rv);
  die_rep_delay_rv = die_rep_delay_rv/sum(die_rep_delay_rv);

  // Cumulative reporting probability
  // for(i in 1:N_days_tot) {
  //   if(i < Max_delay){
  //     cas_cum_report_delay_rv[1+N_days_tot-i] = gamma_cdf(i , cas_rep_delay_shap, cas_rep_delay_rate);
  //     die_cum_report_delay_rv[1+N_days_tot-i] = gamma_cdf(i , die_rep_delay_shap, die_rep_delay_rate);
  //   } else {
  //     cas_cum_report_delay_rv[1+N_days_tot-i] = 1.0;
  //     die_cum_report_delay_rv[1+N_days_tot-i] = 1.0;
  //   }
  // }
  for(i in 1:N_weeks_tot) {
    if(i < Max_delay){
      cas_cum_report_delay_rv[1+N_weeks_tot-i] = gamma_cdf(i , cas_rep_delay_shap, cas_rep_delay_rate);
      die_cum_report_delay_rv[1+N_weeks_tot-i] = gamma_cdf(i , die_rep_delay_shap, die_rep_delay_rate);
    } else {
      cas_cum_report_delay_rv[1+N_weeks_tot-i] = 1.0;
      die_cum_report_delay_rv[1+N_weeks_tot-i] = 1.0;
    }
  }
}
///////////////////////////////////////////////////////////
parameters {
  
// INCIDENCE 
  real                    log_new_inf_0; // starting intercept
  // real<lower=0, upper=6> serial_i; // serial interval
  real serial_i; // serial interval
  // vector[N_spl_par_rt-1]    spl_par_rt0;
  vector[N_spl_par_rt]    spl_par_rt;

// DISEASE PROGRESSION
// probability of transitioning between disease states
  real<lower=0, upper=1>    p_sym_if_inf;
  real<lower=0, upper=1>    p_sev_if_sym;
  real<lower=0, upper=1>    p_die_if_sev;
  real<lower=0>             ifr_decl_OR;
// OMICRON TAKEOVER
// real                        omicron_delay;
  
// DIANGOSIS
// scaling factor for time to diagnosis
  // real<lower=0, upper=1>    scale_dx_delay_asy;
  real<lower=0, upper=1>    scale_dx_delay_sym; 
  real<lower=0, upper=1>    scale_dx_delay_sev; 
// probability of diagnosis at each illness state
  real<lower=0, upper=1>    rr_diag_asy_vs_sym; 
  real<lower=0, upper=1>    p_diag_if_sev;
  vector<lower=0, upper=1>[N_spl_par_dx]  spl_par_sym_dx;

// LIKELIHOOD 
// phi terms for negative b ino imal likelihood function 
  real<lower=0>             inv_sqrt_phi_c;
  real<lower=0>             inv_sqrt_phi_d;
  // VACCINE ADJUSTMENT
  simplex[3]                prob_vac;
}
///////////////////////////////////////////
transformed parameters {
  ///~~~~~~~ Define ~~~~~~~
  // INCIDENCE
  // vector[N_days_tot]      log_new_inf;
  // vector[N_days_tot]      new_inf;
  // vector[N_days_tot]      wane_inf;
  // // vector[N_days_tot]      prot_boost;
  // vector[N_days_tot]      wane_boost;
  // vector[N_days_tot]      deriv1_log_new_inf;
  vector[N_weeks_tot]      log_new_inf;
  vector[N_weeks_tot]      new_inf;
  vector[N_weeks_tot]      wane_inf;
  // vector[N_weeks_tot]      prot_boost;
  vector[N_weeks_tot]      wane_boost;
  vector[N_weeks_tot]      deriv1_log_new_inf;
  real                    pop_uninf;
  real                    pop_sus;
  real                    ever_inf;
  real                    wane_imm;
  // vector[N_spl_par_rt]    spl_par_rt;

  // Rt spline
  // vector[N_days_tot]      logRt0;
  // vector[N_days_tot]      logRt;
  // vector[N_days_tot]      Rt;
  vector[N_weeks_tot]      logRt0;
  vector[N_weeks_tot]      logRt;
  vector[N_weeks_tot]      Rt;
  vector[N_spl_par_rt-1]  deriv1_spl_par_rt;
  vector[N_spl_par_rt-2]  deriv2_spl_par_rt;
  
  // transitions
  vector[N_ifr_adj]   p_die_if_sevt;
  // vector[N_days_tot]  p_sev_if_symt;
  // vector[N_days_tot]  p_sym_if_inft;  
  vector[N_weeks_tot]  p_sev_if_symt;
  vector[N_weeks_tot]  p_sym_if_inft;  
  
  // DIAGNOSIS AND REPORTING  
 // probability of diagnosis
  // vector[N_days_tot]  rr_diag_sym_vs_sev;
  // vector[N_days_tot]  p_diag_if_asy; 
  // vector[N_days_tot]  p_diag_if_sym;
  vector[N_weeks_tot]  rr_diag_sym_vs_sev;
  vector[N_weeks_tot]  p_diag_if_asy; 
  vector[N_weeks_tot]  p_diag_if_sym;

 // daily probabilities of diagnosis and report
 // for days 1 to 60 after entering that state
  vector[Max_delay]  sym_diag_delay_rv;
  vector[Max_delay]  sev_diag_delay_rv;

  // DISEASE OUTCOMES
  // overall case fatality rate
  real                p_die_if_inf;

  // "true" number entering disease state each day
  // vector[N_days_tot]  new_sym; 
  // vector[N_days_tot]  new_sev;
  // vector[N_days_tot]  new_die;
  vector[N_weeks_tot]  new_sym; 
  vector[N_weeks_tot]  new_sev;
  vector[N_weeks_tot]  new_die;

  // newly diagnosed
  // vector[N_days_tot]  new_asy_dx; 
  // vector[N_days_tot]  new_sym_dx; 
  // vector[N_days_tot]  new_sev_dx;
  vector[N_weeks_tot]  new_asy_dx; 
  vector[N_weeks_tot]  new_sym_dx; 
  vector[N_weeks_tot]  new_sev_dx;

  // follow diagnosed cases forward to calculate deaths among diagnosed
  // vector[N_days_tot]  dx_sym_sev; 
  // vector[N_days_tot]  dx_sym_die; 
  // vector[N_days_tot]  dx_sev_die; 
  vector[N_weeks_tot]  dx_sym_sev; 
  vector[N_weeks_tot]  dx_sym_die; 
  vector[N_weeks_tot]  dx_sev_die; 

  // sum to diagnosed cases and deaths
  // vector[N_days_tot]  diag_all;
  // vector[N_days_tot]  new_die_dx;
  vector[N_weeks_tot]  diag_all;
  vector[N_weeks_tot]  new_die_dx;

  // number of cases and deaths in official record on each day
  // (all diagnosed cases with an additional delay to report) 
  // vector[N_days_tot]  occur_cas;
  // vector[N_days_tot]  occur_die; 
  vector[N_weeks_tot]  occur_cas;
  vector[N_weeks_tot]  occur_die; 
  // moving sum cases and deaths
  // vector[N_days_tot]  occur_cas_mvs;
  // vector[N_days_tot]  occur_die_mvs; 
  vector[N_weeks_tot]  occur_cas_mvs;
  vector[N_weeks_tot]  occur_die_mvs; 


  // LIKELIHOOD
  // phi terms for negative binomial likelihood function 
  real phi_cas;
  real phi_die;
 
  // NATURAL HISTORY CASCADE
  p_die_if_sevt = p_die_if_sev * ifr_adj_fixed * (1 + ifr_adj * ifr_decl_OR);

  // for( i in 1:N_days_tot){
  for( i in 1:N_weeks_tot){
  p_die_if_sevt[i]     = p_die_if_sevt[i]   .* pow(ifr_vac_adj[i], prob_vac[1]);
  p_sev_if_symt[i]     = p_sev_if_sym        * pow(ifr_vac_adj[i], prob_vac[2]);
  p_sym_if_inft[i]     = p_sym_if_inf        * pow(ifr_vac_adj[i], prob_vac[3]);
  }
  
  // DIAGNOSIS // 
  // rate ratio of diagnosis at asymptomatic vs symptomatic, symptomatic vs severe
  rr_diag_sym_vs_sev = inv_logit(spl_basis_dx * logit(spl_par_sym_dx));
  
  // probability of diagnosis 
  p_diag_if_sym = p_diag_if_sev * rr_diag_sym_vs_sev;
  p_diag_if_asy = p_diag_if_sym * rr_diag_asy_vs_sym; 
  
  // DELAYS //
  // Diagnosis Delays
  // Calculate the probability of diagnosis for each day in state
  // for 1 to 60 days in state. We do this by scaling the rate term
  // of the gamma distribution of progressiong delays by a modeled fraction 
  // (scale_dx_delay_xxx)
  {
    // Use two vectors to store the results of the `gamma_cdf()` calls in order
    // to avoid double-computing them
    vector[Max_delay+1] sym_delay_gammas;
    vector[Max_delay+1] sev_delay_gammas;
    for (i in 1:Max_delay+1) {
      sym_delay_gammas[i] = gamma_cdf(i-1 , sym_prg_delay_shap, sym_prg_delay_rate/scale_dx_delay_sym);
      sev_delay_gammas[i] = gamma_cdf(i-1 , sev_prg_delay_shap, sev_prg_delay_rate/scale_dx_delay_sev);
    }

    // Diff and reverse, vectorized
    // sym_diag_delay_rv = reverse(sym_delay_gammas[2:(Max_delay+1)] - sym_delay_gammas[1:Max_delay]);
    // sev_diag_delay_rv = reverse(sev_delay_gammas[2:(Max_delay+1)] - sev_delay_gammas[1:Max_delay]);
    for(i in 1:Max_delay){
    sym_diag_delay_rv[1+Max_delay-i] = sym_delay_gammas[i+1] - sym_delay_gammas[i];
    sev_diag_delay_rv[1+Max_delay-i] = sev_delay_gammas[i+1] - sev_delay_gammas[i];
    }

  }

  // DEATHS // 
  // infection fatality rate is the product of the probability of death among
  // severely ill individuals, the probability of being severely ill if 
  // symptomatic, and the probability of becoming symptomatic if infected. 
  p_die_if_inf = p_sym_if_inf * p_sev_if_sym * p_die_if_sev;

  // CASCADE OF INCIDENT OUTCOMES ("TRUE") //

  // NEW INCIDENT CASES
  
  // modeled with a spline
  //     spl_par_rt[2:N_spl_par_rt] = spl_par_rt0; 
  // spl_par_rt[1] = (spl_basis_rt[2+N_days_before,2:N_spl_par_rt]-
  //                  spl_basis_rt[1+N_days_before,2:N_spl_par_rt]) * spl_par_rt0 / 
  //                 (spl_basis_rt[1+N_days_before,1]-spl_basis_rt[2+N_days_before,1]);
  logRt0 = spl_basis_rt * spl_par_rt;

  pop_sus = pop_size * .5;
  pop_uninf = pop_sus;
  ever_inf = 0;

  // for(i in 1:N_days_tot) {
  for(i in 1:N_weeks_tot) {

    logRt[i] = logRt0[i] + log(pop_uninf/pop_size);

    deriv1_log_new_inf[i] = logRt[i]/serial_i;

    log_new_inf[i] = sum(deriv1_log_new_inf[1:i]) + log_new_inf_0;

    new_inf[i] = exp(log_new_inf[i]);

    
    // prot_boost[i] = sum(obs_boost[1:i] * 0.8);
    // if(i > N_days_before){
    // wane_boost[i] = sum(obs_boost[1:i-N_days_before] .* (.35 * exp(-.008 * idx3[N_days_tot-(i-N_days_before)+1:N_days_tot]) + .45) );
    if(i > N_weeks_before){
    wane_boost[i] = sum(obs_boost[1:i-N_weeks_before] .* (.35 * exp(-.008 * idx3[N_weeks_tot-(i-N_weeks_before)+1:N_weeks_tot]*7) + .45) );
    } else{
      wane_boost[i] = obs_boost[1] *.8;
    }
    // wane_imm = .3 * exp(-.008 * idx3[N_days_tot-i+1]) + .2;
    wane_imm = .3 * exp(-.008 * idx3[N_weeks_tot-i+1]*7) + .2;
    pop_sus = pop_size * (1-wane_imm);
    // wane_inf[i] = sum(new_inf[1:i] .* (.75* exp(-.008 * idx3[N_days_tot-i+1:N_days_tot]) + .25));
    wane_inf[i] = sum(new_inf[1:i] .* (.75* exp(-.008 * idx3[N_weeks_tot-i+1:N_weeks_tot] * 7) + .25));
    ever_inf += new_inf[i];
    //CHOOSE ONE OF THE REINFECTION STRATEGIES
    // pop_uninf = pop_sus - ever_inf + prot_boost[i];
    pop_uninf = pop_sus - (wane_inf[i] + wane_boost[i]);
   // pop_uninf -= (new_inf[i] + obs_boost[i]);
   // END OF REINFECTION STRATEGIES
   
    if (pop_uninf < 1) {
      // print("WARNING pop_uninf preliminary value was ", pop_uninf);
      pop_uninf = 1;
    }
  }
  
  Rt = exp(logRt); 
  
  // second derivative
  deriv2_spl_par_rt[1:(N_spl_par_rt-2)] =
    2 * spl_par_rt[2:(N_spl_par_rt-1)] - 
        spl_par_rt[1:(N_spl_par_rt-2)] -
        spl_par_rt[3:N_spl_par_rt]; 
  
  // first derivative
  deriv1_spl_par_rt[1:(N_spl_par_rt-1)] =
    spl_par_rt[2:N_spl_par_rt] - 
    spl_par_rt[1:(N_spl_par_rt-1)];

  // SYMPTOMATIC CASES
  // cases entering a state on day i + j - 1: 
  // cases entering previous state on day i * the probability of progression *
  // the probability progression occurred on day j 
  // print("YO!");
  // print("spl_par_rt:");
  // print(spl_par_rt);
  // print("Rt:");
  // print(Rt);
  // print("new_inf:");
  // print(new_inf);
  // print("inf_prg_delay_rv:");
  // print(inf_prg_delay_rv);
  
  new_sym =
    p_sym_if_inft     .* conv1d(new_inf , inf_prg_delay_rv);

  new_sev = p_sev_if_symt               .* conv1d(new_sym, sym_prg_delay_rv);
  // new_die = p_die_if_sevt[1:N_days_tot] .* conv1d(new_sev, sev_prg_delay_rv);
  new_die = p_die_if_sevt[1:N_weeks_tot] .* conv1d(new_sev, sev_prg_delay_rv);

  // CASCADE OF INCIDENT OUTCOMES (DIAGNOSED) //

  // diagnosed at asymptomatic
  // a diagnosed asymptomatic infection on day i + j - 1
  // is an asymptomatic case on day i with some probability of diagnosis, 
  // and some probability that the diagnosis occurred on day j.
  // we assume asymptomatic diagnosis only occurs among individuals who will be
  // asymptomatic for the entire course of their infection. 
  new_asy_dx = (1 - p_sym_if_inft) .* conv1d(
    new_inf .* p_diag_if_asy,
    asy_rec_delay_rv
  );
  
  // diagnosed at symptomatic
  // a diagnosed symptomatic (not severe) case on day i + j - 1
  // is a symptomatic case on day i with some probability of diagnosis
  // and some probability that the diagnosis occurred on day j.  
  new_sym_dx = conv1d(new_sym .* p_diag_if_sym, sym_diag_delay_rv);
  
  // cascade from diagnosis 
  // follow diagnosed cases forward to determine how many cases diagnosed
  // at symptomatic eventually die 
  dx_sym_sev = p_sev_if_symt .* conv1d(
    new_sym .* p_diag_if_sym,
    sym_prg_delay_rv
  );
        
  // dx_sym_die = p_die_if_sevt[1:N_days_tot] .* conv1d(dx_sym_sev, sev_prg_delay_rv);
  dx_sym_die = p_die_if_sevt[1:N_weeks_tot] .* conv1d(dx_sym_sev, sev_prg_delay_rv);
        
  // diagnosed at severe 
  // as above for symptomatic 
  new_sev_dx = p_diag_if_sev * conv1d(new_sev - dx_sym_sev, sev_diag_delay_rv);
  
  // cascade from diagnosis
  // as above for symptomatic 
  // dx_sev_die = p_diag_if_sev * p_die_if_sevt[1:N_days_tot] .* conv1d(
  //   new_sev - dx_sym_sev,
  //   sev_prg_delay_rv
  // );
  dx_sev_die = p_diag_if_sev * p_die_if_sevt[1:N_weeks_tot] .* conv1d(
    new_sev - dx_sym_sev,
    sev_prg_delay_rv
  );

  // TOTAL DIAGNOSED CASES AND DEATHS //
  diag_all   = new_asy_dx + new_sym_dx + new_sev_dx;
  new_die_dx = dx_sym_die + dx_sev_die;

  // REPORTING //
  // Calcluate "occur_cas" and "occur_die", which are vectors of diagnosed cases 
  // and deaths by the date we expect them to appear in the reported data. 

  // How reporting delays are reflected in the data depend on how the data are 
  // dated.
  //
  // For cases by date of occurrence, we assume all cases diagnosed more than
  // 60 days from the final day of data.
  if(obs_cas_rep == 1)
    occur_cas = conv1d(diag_all, cas_rep_delay_rv);
  else
    occur_cas = diag_all .* cas_cum_report_delay_rv;

  // reporting delays modeled as described above for cases
  if(obs_die_rep == 1)
    occur_die = conv1d(new_sev_dx, die_rep_delay_rv);
  else
    occur_die = new_sev_dx .* die_cum_report_delay_rv;
  // // reporting delays modeled as described above for cases
  // if(obs_die_rep == 1)
  //   occur_die = conv1d(new_die_dx, die_rep_delay_rv);
  // else
  //   occur_die = new_die_dx .* die_cum_report_delay_rv;
    
  // compute moving sums
    // compute the moving sums
  // for(i in 1:N_days_tot) {
  //     if(i < N_days_av) {
  //       occur_cas_mvs[i] = 0;
  //       occur_die_mvs[i] = 0;
  //     } else {
  //       occur_cas_mvs[i] = sum(occur_cas[(i - nda0) : i]);
  //       occur_die_mvs[i] = sum(occur_die[(i - nda0) : i]);
  //     }
  // }
  for(i in 1:N_weeks_tot) {
        occur_cas_mvs[i] = occur_cas[i];
        occur_die_mvs[i] = occur_die[i];
  }

  // phi
  phi_cas = pow(inv_sqrt_phi_c, -2);
  phi_die = pow(inv_sqrt_phi_d, -2);
}
///////////////////////////////////////////////////////////  
model {
  
  // PRIORS
  log_new_inf_0         ~ normal(pri_log_new_inf_0_mu, pri_log_new_inf_0_sd);
  // log_new_inf_0         ~ normal(log_new_inf0calc, pri_log_new_inf_0_sd);
  spl_par_rt            ~ normal(pri_logRt_mu, pri_logRt_sd);
  serial_i              ~ gamma(pri_serial_i_shap, pri_serial_i_rate);
  deriv1_spl_par_rt     ~ normal(0, pri_deriv1_spl_par_sd);
  deriv2_spl_par_rt     ~ normal(0, pri_deriv2_spl_par_sd);

  // PRIORS: DISEASE PROGRESSION
  // probability of transitioning from inf -> sym -> sev -> die
  p_sym_if_inf         ~ beta(pri_p_sym_if_inf_a, pri_p_sym_if_inf_b);
  p_sev_if_sym         ~ beta(pri_p_sev_if_sym_a, pri_p_sev_if_sym_b);
  p_die_if_sev         ~ beta(pri_p_die_if_sev_a, pri_p_die_if_sev_b);
  ifr_decl_OR          ~ gamma(pri_ifr_decl_OR_a, pri_ifr_decl_OR_b);
  
  // PRIORS: overall infection fatality rate
  p_die_if_inf         ~ beta(pri_p_die_if_inf_a, pri_p_die_if_inf_b);

  // PRIORS: diagnosis    
  // probabilities of diagnosis
  rr_diag_asy_vs_sym   ~ beta(pri_rr_diag_asy_vs_sym_a, pri_rr_diag_asy_vs_sym_b);
  spl_par_sym_dx       ~ beta(pri_rr_diag_sym_vs_sev_a,pri_rr_diag_sym_vs_sev_b);
  p_diag_if_sev        ~ beta(pri_p_diag_if_sev_a, pri_p_diag_if_sev_b);

  // delay distribution scaling factors
  scale_dx_delay_sym   ~ beta(scale_dx_delay_sym_a, scale_dx_delay_sym_b); 
  scale_dx_delay_sev   ~ beta(scale_dx_delay_sev_a, scale_dx_delay_sev_b);
  
  // phi  
  inv_sqrt_phi_c       ~ normal(0, 1);
  inv_sqrt_phi_d       ~ normal(0, 1);

  // prop for vaccine
  prob_vac             ~ dirichlet(rep_vector(5, 3));

  // if(N_days_pri_Rt > 0)
  //   logRt[(N_days_tot-N_days_pri_Rt+1) : N_days_tot] ~ normal(0, sd_pri_Rt);
   
  // LIKELIHOOD
  // Before data
  if(pre_period_zero==1){
    // if(N_days_before>0){
    if(N_weeks_before>0){
  
      // if (sum(occur_cas[1:N_days_before]) < 0)
      //   reject("`sum(occur_cas[1:N_days_before])` had a negative value");
      // 
      // if (sum(occur_die[1:N_days_before]) < 0)
      //   reject("`sum(occur_die[1:N_days_before])` had a negative value");
      if (sum(occur_cas[1:N_weeks_before]) < 0)
        reject("`sum(occur_cas[1:N_weeks_before])` had a negative value");

      if (sum(occur_die[1:N_weeks_before]) < 0)
        reject("`sum(occur_die[1:N_weeks_before])` had a negative value");

      // target += neg_binomial_2_lpmf( 0 | sum(occur_cas[1:N_days_before]), phi_cas);
      // target += neg_binomial_2_lpmf( 0 | sum(occur_die[1:N_days_before]), phi_die);
      target += neg_binomial_2_lpmf( 0 | sum(occur_cas[1:N_weeks_before]), phi_cas);
      target += neg_binomial_2_lpmf( 0 | sum(occur_die[1:N_weeks_before]), phi_die);
    }
  }

  if (min(occur_cas) < 0)
    reject("`occur_cas` had a negative value");

  if (min(occur_die) < 0)
    reject("`occur_die` had a negative value");

  // LIKELIHOOD
  // During data
  // target += neg_binomial_2_lpmf(
  //   // `obs_cas` from the first observed day to the last death date
  //   obs_cas_mvs[N_days_av:lastCaseDate] |
  //     // `occur_cas` from the first observed day (`N_days_before+1`) to the
  //     // last death date
  //     occur_cas_mvs[N_days_before+N_days_av : N_days_before+lastCaseDate],
  //   phi_cas
  // ) ;// Optional, but likely unncessesary: / N_days_av;
  // 
  // target += neg_binomial_2_lpmf(
  //   // `obs_die` from the first observed day to the last death date
  //   obs_die[N_days_av:lastDeathDate] |
  //     // `occur_die` from the first observed day (`N_days_before+1`) to the
  //     // last death date
  //     occur_die[N_days_before+N_days_av : N_days_before+lastDeathDate],
  //   phi_die
  // ); // optional, but likelie unnecessary: / N_days_av;
  target += neg_binomial_2_lpmf(
    // `obs_cas` from the first observed day to the last death date
    obs_cas_mvs[1:lastCaseWeek] |
      // `occur_cas` from the first observed day (`N_days_before+1`) to the
      // last death date
      occur_cas_mvs[N_weeks_before+1 : N_weeks_before+lastCaseWeek],
    phi_cas
  ) ;// Optional, but likely unncessesary: / N_days_av;

  target += neg_binomial_2_lpmf(
    // `obs_die` from the first observed day to the last death date
    obs_die_mvs[1:lastDeathWeek] |
      // `occur_die` from the first observed day (`N_days_before+1`) to the
      // last death date
      occur_die_mvs[N_weeks_before+1 : N_weeks_before+lastDeathWeek],
    phi_die
  ); // optional, but likelie unnecessary: / N_days_av;
}
///////////////////////////////////////////////////////////
generated quantities {
  // calculate cumulative incidence + sero_positive + pop_infectiousness
  real                p_die_if_sym;

  // vector[N_days_tot]  diag_cases;  
  vector[N_weeks_tot]  cumulative_incidence;  
  // vector[N_days_tot]  sero_positive;  
  // vector[N_days_tot]  pop_infectiousness;  
  // 
  // vector[Max_delay]   infect_dist_rv;
  // 
  // vector[500]         seropos_dist_rv;

  // cumulative incidence
  cumulative_incidence = cumulative_sum(new_inf); 

  p_die_if_sym = p_die_if_sev * p_sev_if_sym; 

  // diag_cases = new_sym_dx + new_sev_dx; 
  
 // // vector to distribute infectiousness
 //  for(i in 1:Max_delay)
 //    infect_dist_rv[1+Max_delay-i] = 
 //      gamma_cdf(i   , infect_dist_shap, infect_dist_rate) -
 //      gamma_cdf(i-1 , infect_dist_shap, infect_dist_rate);
 //  
 //  // vector to distribute infectiousness seropositives
 //  for(i in 1:500)
 //    seropos_dist_rv[1+500-i] =
 //      1.0 - gamma_cdf(i , seropos_dist_shap, seropos_dist_rate);
 //  
 //  // infectiousness
 //  pop_infectiousness = conv1d(new_inf, infect_dist_rv);
 //  
 //  // seropositives
  // sero_positive = conv1d(new_inf, seropos_dist_rv);
}

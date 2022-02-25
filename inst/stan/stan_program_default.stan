functions {

  // The log beta probability mass of `exp(ltheta)` given alpha parameter `a`
  // and beta parameter `b`.
  //
  // The Stan compiler will automatically generate a function `beta_log_lupdf`
  // which drops the constant additive term `- lbeta(a, b)`. This compiler-
  // -generated function is what will be called in the `model{}` block when
  // incrementing the log-density.
  real beta_log_lpdf(real ltheta, real a, real b) {
    return (a-1)*ltheta + (b-1)*log1m_exp(ltheta) - lbeta(a, b);
  }

  // Included for compatibility with older versions of Stan, which will NOT
  // automatically generate the `_lupdf` version from the `_lpdf` form!
  real beta_log_lupdf(real ltheta, real a, real b) {
    return (a-1)*ltheta + (b-1)*log1m_exp(ltheta);
  }

  vector vlog_sum_exp(vector a, vector b) {
    int l = rows(a);
    vector[l] result;

    if (rows(a) != rows(b))
      reject("nrow(a) must be equal to nrow(b). lengths were: ", rows(a), ", ", rows(b));

    for (i in 1:l)
      result[i] = log_sum_exp(a[i], b[i]);

    return result;
  }

  vector vlog_diff_exp(vector a, vector b) {
    int l = rows(a);
    vector[l] result;

    if (rows(a) != rows(b))
      reject("nrow(a) must be equal to nrow(b). lengths were: ", rows(a), ", ", rows(b));

    for (i in 1:l) {
      if (b[i] >= a[i])
        reject("b[i] was greater than a[i]. b[",i,"]=",b[i],", a[",i,"]=", a[i]);

      result[i] = log_diff_exp(a[i], b[i]);
    }

    return result;
  }

  vector lconv1d(vector lx, vector kernel) {
    int nk = rows(kernel);
    int nx = rows(lx);
    matrix[nx, nk] X; // Matrix-offsets form of the log-space vector `lx`
    matrix[nx, nk] K; // Each row of `K` is the kernel
    matrix[nx, nk] S; // X+K
    vector[nx]     r; // Result vector

    if (nx < nk)
      reject("nrow(x) must be >= nrow(kernel). x had nrow = ", nx);

    // Steps to do a log-space convolution, where `lx` is a vector in
    // log-space, and `kernel` is a vector in real-space:
    //
    // 1. Fill in `X` like usual
    // 2. Create `K` by turning `kernel` into a `row_vector` and using
    //    `rep_matrix()` to make it have the same number of rows as `X`. Each
    //    element of `kernel` is mapped through `log()` so that it's in
    //    log-space.
    // 3. Add `X + K`
    // 4. Call log_sum_exp on each row of the resulting matrix.
    // 5. Return the resulting vector

    // Step 1. Fill in `X` like usual
    //
    // The first `nK-1` rows of the convolution are special, because we don't
    // yet have a full "window" of X entries to take the dot product of against
    // the kernel. We want to handle this case by filling those entires with
    // zeroes. This is NOT to do anything clever; it's to make sure that no
    // entry of the matrix contains undefined values.
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
      X[(nk - i + 1):nx, i] = lx[1:(nx - nk + i)];
    }

    // 2. Create `K` by turning `kernel` into a `row_vector` and using
    //    `rep_matrix()` to make it have the same number of rows as `X`. Each
    //    element of `kernel` is mapped through `log()` so that it's in
    //    log-space.
    K = rep_matrix(log(kernel)', nx);

    // 3. Add `X + K`.
    S = X + K;

    // 4. Call log_sum_exp on each row of the resulting matrix. For the
    //    case where `i < nk`, mask the first `nk - i` entries of the row,
    //    which correspond to the "out-of-range" entries for the convolution.
    for (i in 1:nx) {
      if (i < nk)
        r[i] = log_sum_exp(S[i, (nk - i + 1):nk]);
      else
        r[i] = log_sum_exp(S[i,]);
    }

    // 5. Return the resulting vector
    return r;
  }

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

  real custom_softplus_scalar(real log_val, real log_ceiling) {
    return -1*log_sum_exp(0, -1*(log_val - log_ceiling));
  }

  vector vcustom_softplus_scalar(vector a, vector b) {
    int l = rows(a);
    vector[l] result;

    if (rows(a) != rows(b))
      reject("nrow(a) must be equal to nrow(b). lengths were: ", rows(a), ", ", rows(b));

    for (i in 1:l)
      result[i] = custom_softplus_scalar(a[i], b[i]);

    return result;
  }
}

data {
  // INPUT DATA
  int<lower=0>  N_days;        // days of data
  int<lower=0>  N_days_before; // days before data to init epi model
  int<lower=0>  Max_delay;     // maximum days delay 
  
  int<lower=0>  obs_cas[N_days]; // vector of cases
  int<lower=0>  obs_die[N_days]; // vector of deaths
  real<lower=0> pop_size;        // population size
  
  int<lower=0> N_ifr_adj;                            // length of ifr_adjustment
  vector<lower=0>[N_ifr_adj] ifr_adj;                // ifr_adjustment
  vector<lower=0,upper=1>[N_days+N_days_before] ifr_vac_adj; // ifr_vaccine_adjustment

  real<lower=0> pri_ifr_decl_OR_a; 
  real<lower=0> pri_ifr_decl_OR_b;
  real<lower=0> pri_rr_decl_sev_a;
  real<lower=0> pri_rr_decl_sev_b;
  real<lower=0> pri_rr_decl_die_a;
  real<lower=0> pri_rr_decl_die_b;
  real<lower=0> ifr_adj_fixed;
  
  real<lower=0> infect_dist_rate;
  real<lower=0> infect_dist_shap;
  real<lower=0> seropos_dist_rate;
  real<lower=0> seropos_dist_shap;

  // terms for splines
  // spline parameters and bases
  int<lower=0>                              N_spl_par_rt;
  int<lower=0>                              N_spl_par_dx;
  matrix[N_days+N_days_before,N_spl_par_rt] spl_basis_rt;
  matrix[N_days+N_days_before,N_spl_par_dx] spl_basis_dx;

  // fixed delay distributions
  // time from inf -> sym, sym -> sev, sev -> die
  real<lower=0> inf_prg_delay_shap; 
  real<lower=0> inf_prg_delay_rate;
  real<lower=0> asy_rec_delay_shap;
  real<lower=0> asy_rec_delay_rate;
  real<lower=0> sym_prg_delay_shap;
  real<lower=0> sym_prg_delay_rate;
  real<lower=0> sev_prg_delay_shap;
  real<lower=0> sev_prg_delay_rate;

  // fixed delay from diagnosis to report
  real<lower=0> cas_rep_delay_shap;
  real<lower=0> cas_rep_delay_rate;
  real<lower=0> die_rep_delay_shap;
  real<lower=0> die_rep_delay_rate;
  
  //// control knobs
  // whether to assume zero cases and deaths during warm-up
  int<lower=0, upper=1> pre_period_zero; 

  // what data are included --  cases, deaths: 
  int<lower=0, upper=1> cas_yes; 
  int<lower=0, upper=1> die_yes; 

  // how are the data dated -- report, occurrence: 
  int<lower=0, upper=1> obs_cas_rep;
  int<lower=0, upper=1> obs_die_rep;

  // how many days should be used for the moving average in the likelihood
  // function? 
  int N_days_av; 

  // is there a last obeserved deaths data day?
  int<lower=0> lastDeathDate;

  // is there a last obeserved case data day?
  int<lower=0> lastCaseDate;

  // reinfection setup
  int<lower=0> reinfection;
  int<lower=0> reinf_delay1; 
  int<lower=0> reinf_delay2; 
  vector<lower=0,upper=1>[2] reinf_prob;

  /////////
  // TERMS FOR PRIOR DISTRIBTUIONS
  // for new infections
  real          pri_log_new_inf_0_mu;
  real<lower=0> pri_log_new_inf_0_sd;
  real<lower=0> pri_serial_i_shap; 
  real<lower=0> pri_serial_i_rate; 
  real<lower=0> pri_serial_i_omi_shap; 
  real<lower=0> pri_serial_i_omi_rate; 
  
  // probabilities of progression inf -> sym -> sev -> die
  real<lower=0> pri_p_sym_if_inf_a; 
  real<lower=0> pri_p_sym_if_inf_b;
  real<lower=0> pri_new_p_sym_if_inf_a;
  real<lower=0> pri_new_p_sym_if_inf_b;
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

  // omicron delay
  int<lower=0, upper=1> omicron_adjust;        // 0/1 indicator of whether omicron adjustment should happen
  int                   Omicron_takeover_mean; // ndays from start date that is Dec 20
  real<lower=0>         Omicron_takeover_sd;   // fixed sd for the shape of the omicron takeover. Default: 14
  real<lower=0>         sd_omicron_delay;      // sd of the variation of the mean date: default :10
  
  // input for the number of days to put the Rt prior on
  // int<lower=0>  N_days_pri_Rt;
  // real<lower=0> sd_pri_Rt;
}

///////////////////////////////////////////////////////////
transformed data {
  int N_days_tot;

  // Moving sums
  int<lower=0> obs_cas_mvs[N_days]; // vector of cases
  int<lower=0> obs_die_mvs[N_days]; // vector of deaths
  int<lower=0> nda0 = N_days_av - 1;

  vector<upper=0>[N_days+N_days_before] log_ifr_vac_adj; // ifr_vaccine_adjustment

  real log_ifr_adj_fixed;
  real log_pop_size;

  // Kernel for doing 7-day sums of quantities.
  vector<lower=0>[nda0] sum_nda0_kernel;

  // Progression delays
  vector<lower=0>[Max_delay] inf_prg_delay_rv;
  vector<lower=0>[Max_delay] asy_rec_delay_rv; 
  vector<lower=0>[Max_delay] sym_prg_delay_rv;
  vector<lower=0>[Max_delay] sev_prg_delay_rv;
 
  // Reporting delays
  vector<lower=0>[Max_delay] cas_rep_delay_rv;
  vector<lower=0>[Max_delay] die_rep_delay_rv;
 
  // Cumulative reporting delays
  vector<lower=0>[N_days + N_days_before] cas_cum_report_delay_rv; 
  vector<lower=0>[N_days + N_days_before] die_cum_report_delay_rv; 
  vector[N_days + N_days_before] log_cas_cum_report_delay_rv; 
  vector[N_days + N_days_before] log_die_cum_report_delay_rv; 
  
  // create 'N_days_tot', which is days of data plus days to model before first 
  // case or death 
  N_days_tot = N_days + N_days_before; 

  // compute the moving sums
  for(i in 1:N_days) {
    if(i < N_days_av) {
      obs_cas_mvs[i] = 0;
      obs_die_mvs[i] = 0;
    } else {
      obs_cas_mvs[i] = sum(obs_cas[(i - nda0) : i]);
      obs_die_mvs[i] = sum(obs_die[(i - nda0) : i]);
    }
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

  log_ifr_vac_adj   = log(ifr_vac_adj);
  log_ifr_adj_fixed = log(ifr_adj_fixed);
  log_pop_size      = log(pop_size);
  
  sum_nda0_kernel = rep_vector(1, nda0);

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
  for(i in 1:N_days_tot) {
    if(i < Max_delay){
      cas_cum_report_delay_rv[1+N_days_tot-i] = gamma_cdf(i, cas_rep_delay_shap, cas_rep_delay_rate);
      die_cum_report_delay_rv[1+N_days_tot-i] = gamma_cdf(i, die_rep_delay_shap, die_rep_delay_rate);
    } else {
      cas_cum_report_delay_rv[1+N_days_tot-i] = 1.0;
      die_cum_report_delay_rv[1+N_days_tot-i] = 1.0;
    }
  }

  log_cas_cum_report_delay_rv = log(cas_cum_report_delay_rv);
  log_die_cum_report_delay_rv = log(die_cum_report_delay_rv);
}
///////////////////////////////////////////////////////////
parameters {
  
  // INCIDENCE 
  // Bounds are for BFGS.
  real<lower=-15, upper=log(1000)> log_new_inf_0; // starting intercept
  real<lower=1, upper=11>   serial_i; // serial interval
  // real<lower=1, upper=8>    serial_i_omi; // serial interval

  // The bounds on this vector keep Rt strictly positive and also prevent
  // any completely insane Rt values from occurring. This seems to help
  // greatly with initialization
  vector<lower=0.001, upper=5>[N_spl_par_rt] spl_par_rt_raw;

  // DISEASE PROGRESSION
  // probability of transitioning between disease states
  real<upper=0> log_p_sym_if_inf;
  real<upper=0> log_p_sym_if_inf_omi;
  real<upper=0> log_p_sev_if_sym;
  real<upper=0> log_p_die_if_sev;

  real<lower=0.001, upper=4> ifr_decl_OR;
  // Bounds are for BFGS.
  real<lower=0.01, upper=0.3> rr_decl_sev;
  // Bounds are for BFGS.
  real<lower=0.0001, upper=0.25> rr_decl_die;

  // OMICRON TAKEOVER
  real<lower=-60, upper=60> omicron_delay;
  
  // DIANGOSIS
  // scaling factor for time to diagnosis
  real<upper=0> log_scale_dx_delay_sym; 
  real<upper=0> log_scale_dx_delay_sev; 

  // probability of diagnosis at each illness state
  real<upper=0> log_rr_diag_asy_vs_sym; 
  real<upper=0> log_p_diag_if_sev;

  vector<lower=0.001, upper=0.999>[N_spl_par_dx] spl_par_sym_dx;

  // LIKELIHOOD 
  // phi terms for negative binomimal likelihood function 
  real<lower=0.001, upper=15> inv_sqrt_phi_c;
  real<lower=0.001, upper=15> inv_sqrt_phi_d;

  // reinfection probability 
  // Bounds are for BFGS.
  real<upper=0> log_p_reinf;

  // VACCINE ADJUSTMENT
  // Bounds are for BFGS.
  simplex[3] prob_vac;
}

///////////////////////////////////////////
transformed parameters {
  ///~~~~~~~ Define ~~~~~~~
  // INCIDENCE
  vector[N_days_tot] log_new_inf;

  // The log of the fraction of the population that is susceptible.
  real<upper=0> log_frac_sus;

  // The fraction of the susceptible population that becomes infected at
  // a given timestep. This is only used in the infection loop.
  //
  // The denominator is effectively `exp(log_frac_sus)`.
  real<upper=0> log_frac_new_inf;

  // The number of susceptible individuals in the population at any given
  // point in the loop.
  real log_pop_sus;

  vector[N_days_tot] log_serial_i_comb;

  // Rt spline - the product of the spline design matrix and the spline
  // parameters.
  vector<lower=0>[N_days_tot] Rt0;
  vector[N_days_tot] logRt0;
  vector[N_days_tot] logRt;
  
  // Transitions.
  //
  // Whether or not some more of these can be upper-bounded at 1 should be
  // reviewed.
  vector<upper=log(4)>[N_days_tot] log_p_die_if_sevt;
  vector<upper=log(4)>[N_days_tot] log_p_sev_if_symt;
  vector<upper=0>[N_days_tot] log_p_sym_if_inft;  
  vector<upper=0>[N_days_tot] log_p_sym_if_inft_omi;  

  // new probability of symptomatic
  real log_rr_sym_if_inf;
  real log_rr_sev_if_sym;
  real<lower=0,upper=1> rr_die_if_sev;
  
  // DIAGNOSIS AND REPORTING  
  // probability of diagnosis
  vector[N_days_tot] rr_diag_sym_vs_sev;
  vector<upper=0>[N_days_tot] log_p_diag_if_asy; 
  vector<upper=0>[N_days_tot] log_p_diag_if_sym;

  // daily probabilities of diagnosis and report
  // for days 1 to 60 after entering that state
  vector<lower=0>[Max_delay] sym_diag_delay_rv;
  vector<lower=0>[Max_delay] sev_diag_delay_rv;

  // DISEASE OUTCOMES
  // overall case fatality rate
  real<upper=0> log_p_die_if_inf;

  // "true" number entering disease state each day
  // vector[N_days_tot] new_sym; 
  vector[N_days_tot] log_new_sym; 
  // vector[N_days_tot] new_sev;
  vector[N_days_tot] log_new_sev;
  // vector[N_days_tot] new_die;
  vector[N_days_tot] log_new_die;

  // newly diagnosed
  // vector[N_days_tot] new_asy_dx; 
  vector[N_days_tot] log_new_asy_dx; 
  // vector[N_days_tot] new_sym_dx; 
  vector[N_days_tot] log_new_sym_dx; 
  // vector[N_days_tot] new_sev_dx;
  vector[N_days_tot] log_new_sev_dx;

  // follow diagnosed cases forward to calculate deaths among diagnosed
  // vector[N_days_tot] dx_sym_sev; 
  vector[N_days_tot] log_dx_sym_sev; 
  // vector[N_days_tot] dx_sym_die; 
  vector[N_days_tot] log_dx_sym_die; 
  // vector[N_days_tot] dx_sev_die; 
  vector[N_days_tot] log_dx_sev_die; 

  // sum to diagnosed cases and deaths
  // vector[N_days_tot] diag_all;
  vector[N_days_tot] log_diag_all;
  // vector[N_days_tot] new_die_dx;
  vector[N_days_tot] log_new_die_dx;

  // number of cases and deaths in official record on each day
  // (all diagnosed cases with an additional delay to report) 
  // vector[N_days_tot] occur_cas;
  vector[N_days_tot] log_occur_cas;
  // vector[N_days_tot] occur_die; 
  vector[N_days_tot] log_occur_die; 

  // moving sum cases and deaths
  // vector[N_days_tot] occur_cas_mvs;
  // vector[N_days_tot] occur_die_mvs; 
  vector[N_days_tot] log_occur_cas_mvs;
  vector[N_days_tot] log_occur_die_mvs; 

  // OMICRON DELAY int 
  vector<upper=0>[N_days_tot] log_ifr_omi_rv;
  vector<upper=0>[N_days_tot] log_ifr_omi_rv_sev;
  vector<upper=0>[N_days_tot] log_ifr_omi_rv_die;

  // LIKELIHOOD
  // phi terms for negative binomial likelihood function 
  real<lower=0> phi_cas;
  real<lower=0> phi_die;

  // ========== PRINT PARAMETERS ==============
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  //
  // print("log_new_inf_0:");
  // print(log_new_inf_0);
  // print("serial_i:");
  // print(serial_i);
  // print("serial_i_omi:");
  // print(serial_i_omi);
  // print("spl_par_rt_raw:");
  // print(spl_par_rt_raw);
  // print("p_sym_if_inf:");
  // print(p_sym_if_inf);
  // print("p_sym_if_inf_omi:");
  // print(p_sym_if_inf_omi);
  // print("p_sev_if_sym:");
  // print(p_sev_if_sym);
  // print("p_die_if_sev:");
  // print(p_die_if_sev);
  // print("ifr_decl_OR:");
  // print(ifr_decl_OR);
  // print("rr_decl_sev:");
  // print(rr_decl_sev);
  // print("rr_decl_die:");
  // print(rr_decl_die);
  // print("omicron_delay:");
  // print(omicron_delay);
  // print("scale_dx_delay_sym:");
  // print(scale_dx_delay_sym);
  // print("scale_dx_delay_sev:");
  // print(scale_dx_delay_sev);
  // print("rr_diag_asy_vs_sym:");
  // print(rr_diag_asy_vs_sym);
  // print("p_diag_if_sev:");
  // print(p_diag_if_sev);
  // print("spl_par_sym_dx:");
  // print(spl_par_sym_dx);
  // print("inv_sqrt_phi_c:");
  // print(inv_sqrt_phi_c);
  // print("inv_sqrt_phi_d:");
  // print(inv_sqrt_phi_d);
  // print("p_reinf:");
  // print(p_reinf);
  // print("prob_vac:");
  // print(prob_vac);
  //
  //
  //
  //
  //
  //
  //
  //
  //
  // ====== END PRINT PARAMETERS ==============

  // Omicron-takeover sigmoid. First defined in realspace, then transformed to
  // logspace.
  //
  // For small values of `i`, `normal_cdf` underflows to 0. Additionally, we
  // believe that Omicron never takes over 100%. So, we scale the vectors by
  // `0.95`, then we translate them by 0.001. This ensures that the domain of
  // the vectors is [0.001, 0.951]. This guarantees that when we use `log()` on
  // the vectors, they do not go to `-Inf`.
  log_ifr_omi_rv     = rep_vector(0, N_days_tot);
  log_ifr_omi_rv_die = rep_vector(0, N_days_tot);
  log_ifr_omi_rv_sev = rep_vector(0, N_days_tot);

  for (i in (Omicron_takeover_mean - 50):N_days_tot) {
    log_ifr_omi_rv[i]     = normal_cdf(i, Omicron_takeover_mean + omicron_delay,             Omicron_takeover_sd);
    log_ifr_omi_rv_die[i] = normal_cdf(i, Omicron_takeover_mean + omicron_delay + 6 + 7 + 9, Omicron_takeover_sd);
    log_ifr_omi_rv_sev[i] = normal_cdf(i, Omicron_takeover_mean + omicron_delay + 6 + 7,     Omicron_takeover_sd);
  }

  log_ifr_omi_rv     *= 0.95;
  log_ifr_omi_rv_die *= 0.95;
  log_ifr_omi_rv_sev *= 0.95;
  log_ifr_omi_rv     += 0.001;
  log_ifr_omi_rv_die += 0.001;
  log_ifr_omi_rv_sev += 0.001;
  log_ifr_omi_rv     = log(log_ifr_omi_rv);
  log_ifr_omi_rv_die = log(log_ifr_omi_rv_die);
  log_ifr_omi_rv_sev = log(log_ifr_omi_rv_sev);
  
  // compute new serial i, combination of ancestral and omicron
  log_serial_i_comb[1:(Omicron_takeover_mean - 50 - 1)] =
    rep_vector(log(serial_i), Omicron_takeover_mean  - 50 - 1);

  log_serial_i_comb[(Omicron_takeover_mean - 50):N_days_tot] =
    vlog_sum_exp(
      log(serial_i) +
        log1m_exp(log_ifr_omi_rv[(Omicron_takeover_mean - 50):N_days_tot]),
      log(2.5) +
        log_ifr_omi_rv[(Omicron_takeover_mean - 50):N_days_tot] // Fixed omicron serial interval @ 2.5 days!
    );
  
  // RELATIVE RISKS for omicron adjustment 
  log_rr_sym_if_inf = log_p_sym_if_inf_omi - log_p_sym_if_inf;
  log_rr_sev_if_sym = log(rr_decl_sev) - log_rr_sym_if_inf;
  rr_die_if_sev = rr_decl_die / rr_decl_sev;

  // NATURAL HISTORY CASCADE
  // p_die_if_sevt = p_die_if_sev * ifr_adj_fixed * (1 + ifr_adj * ifr_decl_OR);
  log_p_die_if_sevt = log_p_die_if_sev + log_ifr_adj_fixed +
    log1p(ifr_adj[:N_days_tot] * ifr_decl_OR);

  // log_p_die_if_sevt[i] = p_die_if_sevt[i] *
  //   pow(ifr_vac_adj[i], prob_vac[1]) *
  //   (1 - exp(log_ifr_omi_rv_die[i]) * (1 - rr_die_if_sev));
  log_p_die_if_sevt = log_p_die_if_sevt +
    prob_vac[1] * log_ifr_vac_adj +
    log1m_exp(log_ifr_omi_rv_die + log1m(rr_die_if_sev));

  // p_sev_if_symt[i] = p_sev_if_sym *
  //   pow(ifr_vac_adj[i], prob_vac[2]) *
  //   (1 - exp(log_ifr_omi_rv_sev[i]) * (1 - rr_sev_if_sym));
  log_p_sev_if_symt = log_p_sev_if_sym +
    prob_vac[2] * log_ifr_vac_adj +
    log1m_exp(log_ifr_omi_rv_sev + log1m_exp(log_rr_sev_if_sym));

  log_p_sym_if_inft     = log_p_sym_if_inf     + prob_vac[3]*log_ifr_vac_adj;
  log_p_sym_if_inft_omi = log_p_sym_if_inf_omi + prob_vac[3]*log_ifr_vac_adj;
  
  // DIAGNOSIS // 
  // rate ratio of diagnosis at asymptomatic vs symptomatic, symptomatic vs severe
  rr_diag_sym_vs_sev = inv_logit(spl_basis_dx * logit(spl_par_sym_dx));
  
  // probability of diagnosis 
  log_p_diag_if_sym = log_p_diag_if_sev + log(rr_diag_sym_vs_sev);
  log_p_diag_if_asy = log_p_diag_if_sym + log_rr_diag_asy_vs_sym; 
  
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
      sym_delay_gammas[i] = gamma_cdf(i-1, sym_prg_delay_shap, sym_prg_delay_rate/exp(log_scale_dx_delay_sym));
      sev_delay_gammas[i] = gamma_cdf(i-1, sev_prg_delay_shap, sev_prg_delay_rate/exp(log_scale_dx_delay_sev));
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
  log_p_die_if_inf = log_p_sym_if_inf + log_p_sev_if_sym + log_p_die_if_sev;

  // CASCADE OF INCIDENT OUTCOMES ("TRUE") //

  // NEW INCIDENT CASES
  
  // modeled with a spline
  Rt0 = spl_basis_rt * spl_par_rt_raw;
  logRt0 = log(Rt0);
  log_pop_sus = log_pop_size;

  // `log_frac_sus` can't be defined as 1, because it's passed directly
  // to `log1m_exp()` during the course of the infection loop. Instead we
  // choose a sufficiently similar number so that with the typical population
  // size of geographies, we will be off by fewer than 1 person.
  log_frac_sus = log(0.99999999999);

  for(i in 1:N_days_tot) {
    logRt[i] = logRt0[i] + log_pop_sus - log_pop_size;

    if (i > 1) {
      // ** NON-LOGSPACE IMPL **
      // new_inf[i] = new_inf[i-1] * Rt[i]^(1/serial_i_comb[i]);
      // ** LOGSPACE IMPL **
      log_frac_new_inf = custom_softplus_scalar(
        log_new_inf[i-1] + (1/exp(log_serial_i_comb[i])) * logRt[i],
        log_pop_sus
      );

      log_new_inf[i] = log_frac_new_inf + log_pop_sus;
    }
    else if (i == 1) {
      // ** NON-LOGSPACE IMPL **
      // new_inf[i] = Rt[i]^(1/serial_i_comb[i]) * exp(log_new_inf_0);
      //
      // ** LOGSPACE IMPL **
      log_frac_new_inf = custom_softplus_scalar(

        // This could underflow and lead to a divide-by-0 problem, however,
        // we have not yet seen this occur in practice
        log_new_inf_0 + (1/exp(log_serial_i_comb[i])) * logRt[i],
        log_pop_sus
      );

      log_new_inf[i] = log_frac_new_inf + log_pop_sus;
    }
    
    // The fraction of the population that's susceptible
    log_frac_sus = log_sum_exp(

      // This first term is the fraction of the population that
      // remains after subtracting out the susceptible individuals who were
      // just infected.
      //
      // Stated as: one minus the fraction of the susceptibles that were just
      // infected, times the old log_frac_sus.
      //
      // This "works" because, as a fraction, we can constrain
      // `log_frac_new_inf` to be < 0. The implicit "denominator" for
      // `log_frac_new_inf` is the "old" `log_frac_sus`.
      log_frac_sus + log1m_exp(log_frac_new_inf),
    
      // Second term is the fraction of the "non-susceptible" population that
      // gets moved to the "susceptible" category.
      log_p_reinf + log_ifr_omi_rv[i] + log1m_exp(log_frac_sus)
    );

    log_pop_sus = log_frac_sus + log_pop_size;
  }

  // SYMPTOMATIC CASES
  // cases entering a state on day i + j - 1: 
  // cases entering previous state on day i * the probability of progression *
  // the probability progression occurred on day j 

  // ** NON-LOGSPACE IMPL **
  // new_sym =
  //   p_sym_if_inft     .* conv1d(new_inf .* (1 - ifr_omi_rv), inf_prg_delay_rv) +
  //   p_sym_if_inft_omi .* conv1d(new_inf .* ifr_omi_rv      , inf_prg_delay_rv);
  //
  // ** LOGSPACE IMPL **
  log_new_sym = vlog_sum_exp(
    log_p_sym_if_inft     + lconv1d(log_new_inf + log1m_exp(log_ifr_omi_rv), inf_prg_delay_rv),
    log_p_sym_if_inft_omi + lconv1d(log_new_inf + log_ifr_omi_rv,            inf_prg_delay_rv)
  );

  // ** NON-LOGSPACE IMPL **
  // new_sev = p_sev_if_symt .* conv1d(new_sym, sym_prg_delay_rv);
  //
  // ** LOGSPACE IMPL **
  log_new_sev = log_p_sev_if_symt + lconv1d(log_new_sym, sym_prg_delay_rv);


  // ** NON-LOGSPACE IMPL **
  // new_die = p_die_if_sevt[1:N_days_tot] .* conv1d(new_sev, sev_prg_delay_rv);
  //
  // ** LOGSPACE IMPL **
  log_new_die = log_p_die_if_sevt + lconv1d(log_new_sev, sev_prg_delay_rv);

  // CASCADE OF INCIDENT OUTCOMES (DIAGNOSED) //

  // Diagnosed at asymptomatic
  //
  // A diagnosed asymptomatic infection on day i + j - 1 is an asymptomatic
  // case on day i with some probability of diagnosis, and some probability
  // that the diagnosis occurred on day j.  we assume asymptomatic diagnosis
  // only occurs among individuals who will be asymptomatic for the entire
  // course of their infection. 
  //
  // ** NON-LOGSPACE IMPL **
  // new_asy_dx = (1 - p_sym_if_inft) .* conv1d(
  //   new_inf .* p_diag_if_asy,
  //   asy_rec_delay_rv
  // );
  //
  // ** LOGSPACE IMPL **
  log_new_asy_dx = log1m_exp(log_p_sym_if_inft) + lconv1d(
    log_new_inf + log_p_diag_if_asy,
    asy_rec_delay_rv
  );
  
  // Diagnosed at symptomatic
  //
  // A diagnosed symptomatic (not severe) case on day i + j - 1
  // is a symptomatic case on day i with some probability of diagnosis
  // and some probability that the diagnosis occurred on day j.  
  //
  // ** NON-LOGSPACE IMPL **
  // new_sym_dx = conv1d(new_sym .* p_diag_if_sym, sym_diag_delay_rv);
  //
  // ** LOGSPACE IMPL **
  log_new_sym_dx = lconv1d(log_p_diag_if_sym + log_new_sym, sym_diag_delay_rv);
  
  // Cascade from diagnosis 
  //
  // Follow diagnosed cases forward to determine how many cases diagnosed
  // at symptomatic eventually die.
  //
  // ** NON-LOGSPACE IMPL **
  // dx_sym_sev = p_sev_if_symt .* conv1d(
  //   new_sym .* p_diag_if_sym,
  //   sym_prg_delay_rv
  // );
  //
  // ** LOGSPACE IMPL **
  log_dx_sym_sev = log_p_sev_if_symt + lconv1d(
    log_p_diag_if_sym + log_new_sym,
    sym_prg_delay_rv
  );
  
  // ** NON-LOGSPACE IMPL **
  // dx_sym_die = p_die_if_sevt[1:N_days_tot] .* conv1d(dx_sym_sev, sev_prg_delay_rv);
  //
  // ** LOGSPACE IMPL **
  log_dx_sym_die =
    log_p_die_if_sevt + lconv1d(log_dx_sym_sev, sev_prg_delay_rv);
        
  // Diagnosed at severe 
  //
  // As above for symptomatic 
  //
  // ** NON-LOGSPACE IMPL **
  // new_sev_dx = p_diag_if_sev * conv1d(new_sev - dx_sym_sev, sev_diag_delay_rv);
  //
  // ** LOGSPACE IMPL **
  log_new_sev_dx = log_p_diag_if_sev +
    lconv1d(
      vlog_diff_exp(
        log_new_sev,

        // NEW: Clamping function, to avoid underflow condition. May not work
        // as intended! Further work needed to verify that expected values of
        // `log_dx_sym_sev` never reach the region of this function where the
        // clamping behavior becomes strong.
        log_new_sev + vcustom_softplus_scalar(log_dx_sym_sev, log_new_sev)
      ),
      sev_diag_delay_rv
    );
  
  // Cascade from diagnosis
  //
  // As above for symptomatic.
  //
  // ** NON-LOGSPACE IMPL **
  // dx_sev_die = p_die_if_sev * p_die_if_sevt[1:N_days_tot] .* conv1d(
  //   new_sev - dx_sym_sev,
  //   sev_prg_delay_rv
  // );
  //
  // ** LOGSPACE IMPL **
  log_dx_sev_die =
    log_p_die_if_sev +
    log_p_die_if_sevt +
    lconv1d(
      vlog_diff_exp(
        log_new_sev,

        // NEW: Clamping function, to avoid underflow condition. May not work
        // as intended! Further work needed to verify that expected values of
        // `log_dx_sym_sev` never reach the region of this function where the
        // clamping behavior becomes strong.
        log_new_sev + vcustom_softplus_scalar(log_dx_sym_sev, log_new_sev)
      ),
      sev_prg_delay_rv
    );

  // TOTAL DIAGNOSED CASES AND DEATHS //
  //
  // ** NON-LOGSPACE IMPL **
  // diag_all   = new_asy_dx + new_sym_dx + new_sev_dx;
  //
  // ** LOGSPACE IMPL **
  log_diag_all = vlog_sum_exp(
    log_new_asy_dx,
    vlog_sum_exp(
      log_new_sym_dx,
      log_new_sev_dx
    )
  );

  // ** NON-LOGSPACE IMPL **
  // new_die_dx = dx_sym_die + dx_sev_die;
  //
  // ** LOGSPACE IMPL **
  log_new_die_dx = vlog_sum_exp(log_dx_sym_die, log_dx_sev_die);

  // REPORTING //
  // Calcluate "occur_cas" and "occur_die", which are vectors of diagnosed cases 
  // and deaths by the date we expect them to appear in the reported data. 

  // How reporting delays are reflected in the data depend on how the data are 
  // dated.
  //
  // For cases by date of occurrence, we assume all cases diagnosed more than
  // 60 days from the final day of data.
  if(obs_cas_rep == 1) {                              
    // ** NON-LOGSPACE IMPL **
    // occur_cas = conv1d(diag_all, cas_rep_delay_rv);
    //
    // ** LOGSPACE IMPL **
    log_occur_cas = lconv1d(log_diag_all, cas_rep_delay_rv);
  } else {
    // ** NON-LOGSPACE IMPL **
    // occur_cas = diag_all .* cas_cum_report_delay_rv;
    //
    // ** LOGSPACE IMPL **
    log_occur_cas = log_cas_cum_report_delay_rv + log_diag_all;
  }

  // reporting delays modeled as described above for cases
  if(obs_die_rep == 1) {
    // ** NON-LOGSPACE IMPL **
    // occur_die = conv1d(new_die_dx, die_rep_delay_rv);
    //
    // ** LOGSPACE IMPL **
    log_occur_die = lconv1d(log_new_die_dx, die_rep_delay_rv);
  } else {
    // ** NON-LOGSPACE IMPL **
    // occur_die = new_die_dx .* die_cum_report_delay_rv;
    //
    // ** LOGSPACE IMPL **
    log_occur_die = log_die_cum_report_delay_rv + log_new_die_dx;
  }
    
  // Convolve the occur* vectors with (7) 1's - this is equivalent to summing
  // the last `nda0` (7) days of entries.
  log_occur_cas_mvs = lconv1d(log_occur_cas, sum_nda0_kernel);
  log_occur_die_mvs = lconv1d(log_occur_die, sum_nda0_kernel);

  // phi
  phi_cas = pow(inv_sqrt_phi_c, -2);
  phi_die = pow(inv_sqrt_phi_d, -2);
}
///////////////////////////////////////////////////////////  
model {
  
  // PRIORS
  log_new_inf_0        ~ normal(pri_log_new_inf_0_mu, pri_log_new_inf_0_sd);
                                   
  spl_par_rt_raw       ~ gamma(4.4, 2.9);

  spl_par_rt_raw[2:] - spl_par_rt_raw[:N_spl_par_rt-1] ~ std_normal();

  serial_i             ~ gamma(pri_serial_i_shap, pri_serial_i_rate);

  // Currently FIXING the Omicron serial interval for testing
  // serial_i_omi         ~ gamma(pri_serial_i_omi_shap, pri_serial_i_omi_rate);

  // PRIORS: DISEASE PROGRESSION
  //
  // Probability of transitioning from inf -> sym -> sev -> die
  target += beta_log_lupdf(log_p_sym_if_inf, pri_p_sym_if_inf_a, pri_p_sym_if_inf_b);

  // Currently TESTING a tighter prior on `p_sym_if_inf_omi`.
  //
  // p_sym_if_inf_omi     ~ beta(pri_new_p_sym_if_inf_a, pri_new_p_sym_if_inf_b);
  target += beta_log_lupdf(log_p_sym_if_inf_omi, 5, 50);
  target += beta_log_lupdf(log_p_sev_if_sym, pri_p_sev_if_sym_a, pri_p_sev_if_sym_b);
  target += beta_log_lupdf(log_p_die_if_sev, pri_p_die_if_sev_a, pri_p_die_if_sev_b);

  // PRIORS: Overall infection fatality rate
  target += beta_log_lupdf(log_p_die_if_inf, pri_p_die_if_inf_a, pri_p_die_if_inf_b);

  // PRIORS: Diagnosis    
  //
  // Probabilities of diagnosis
  target += beta_log_lupdf(log_rr_diag_asy_vs_sym, pri_rr_diag_asy_vs_sym_a, pri_rr_diag_asy_vs_sym_b);
  spl_par_sym_dx ~ beta(pri_rr_diag_sym_vs_sev_a, pri_rr_diag_sym_vs_sev_b);
  target += beta_log_lupdf(log_p_diag_if_sev, pri_p_diag_if_sev_a, pri_p_diag_if_sev_b);

  // Delay distribution scaling factors
  target += beta_log_lupdf(log_scale_dx_delay_sym, scale_dx_delay_sym_a, scale_dx_delay_sym_b); 
  target += beta_log_lupdf(log_scale_dx_delay_sev, scale_dx_delay_sev_a, scale_dx_delay_sev_b);
  
  ifr_decl_OR          ~ gamma(pri_ifr_decl_OR_a, pri_ifr_decl_OR_b);
  rr_decl_sev          ~ gamma(pri_rr_decl_sev_a, pri_rr_decl_sev_b);
  rr_decl_die          ~ gamma(pri_rr_decl_die_a, pri_rr_decl_die_b);
  
  // Phi  
  inv_sqrt_phi_c       ~ std_normal();
  inv_sqrt_phi_d       ~ std_normal();
  
  // Reinfection probability
  target += beta_log_lupdf(log_p_reinf, 4,32);

  // Omicron delay 
  omicron_delay        ~ normal(0, sd_omicron_delay);

  // Prop for vaccine
  prob_vac             ~ dirichlet(rep_vector(5, 3));

  // LIKELIHOOD
  //
  // Before data
  
  // ** NON-LOGSPACE IMPL **
  // target += neg_binomial_2_lpmf( 0 | sum(occur_cas[1:N_days_before]), phi_cas);
  // target += neg_binomial_2_lpmf( 0 | sum(occur_die[1:N_days_before]), phi_die);
  //
  // ** LOGSPACE IMPL **
  if(N_days_av < N_days_before){
    target += neg_binomial_2_log_lpmf(
      0 | log_sum_exp(log_occur_cas[N_days_av:N_days_before]),
      phi_cas
    );

    target += neg_binomial_2_log_lpmf(
      0 | log_sum_exp(log_occur_die[N_days_av:N_days_before]),
      phi_die
    );
  }

  // LIKELIHOOD
  // During data

  // ** NON-LOGSPACE IMPL **
  // target += neg_binomial_2_lpmf(
  //   // `obs_cas` from the first observed day to the last death date
  //   obs_cas_mvs[N_days_av:lastCaseDate] |
  //     // `occur_cas` from the first observed day (`N_days_before+1`) to the
  //     // last death date
  //     occur_cas_mvs[N_days_before+N_days_av : N_days_before+lastCaseDate],
  //   phi_cas
  // ) ;// Optional, but likely unncessesary: / N_days_av;
  //
  // ** LOGSPACE IMPL **
  target += neg_binomial_2_log_lpmf(
    // `obs_cas` from the first observed day to the last death date
    obs_cas_mvs[N_days_av:lastCaseDate] |
      // `occur_cas` from the first observed day (`N_days_before+1`) to the
      // last death date
      log_occur_cas_mvs[N_days_before+N_days_av : N_days_before+lastCaseDate],
    phi_cas
  ) ;// Optional, but likely unncessesary: `- log(N_days_av)`;

  // ** NON-LOGSPACE IMPL **
  // target += neg_binomial_2_log_lpmf(
  //   // `obs_die` from the first observed day to the last death date
  //   obs_die[N_days_av:lastDeathDate] |
  //     // `occur_die` from the first observed day (`N_days_before+1`) to the
  //     // last death date
  //     occur_die[N_days_before+N_days_av : N_days_before+lastDeathDate],
  //   phi_die
  // ); // optional, but likely unnecessary: / N_days_av;
  //
  // ** LOGSPACE IMPL **
  target += neg_binomial_2_log_lpmf(
    // `obs_die` from the first observed day to the last death date
    obs_die[N_days_av:lastDeathDate] |
      // `occur_die` from the first observed day (`N_days_before+1`) to the
      // last death date
      log_occur_die[N_days_before+N_days_av : N_days_before+lastDeathDate],
    phi_die
  ); // optional, but likely unnecessary: / N_days_av;
}
///////////////////////////////////////////////////////////
generated quantities {
  // calculate cumulative incidence + sero_positive + pop_infectiousness
  // real                p_die_if_sym;

  // vector[N_days_tot]  diag_cases;  
  // vector[N_days_tot]  cumulative_incidence;  
  // vector[N_days_tot]  sero_positive;  
  // vector[N_days_tot]  pop_infectiousness;  

  // vector[Max_delay]   infect_dist_rv;

  // vector[500]         seropos_dist_rv;

  // cumulative incidence
  // cumulative_incidence = cumulative_sum(new_inf); 

  // p_die_if_sym = p_die_if_sev * p_sev_if_sym; 

  // diag_cases = new_sym_dx + new_sev_dx; 
  // 
  // // vector to distribute infectiousness
  // for(i in 1:Max_delay)
  //   infect_dist_rv[1+Max_delay-i] = 
  //     gamma_cdf(i   , infect_dist_shap, infect_dist_rate) -
  //     gamma_cdf(i-1 , infect_dist_shap, infect_dist_rate);
  // 
  // // vector to distribute infectiousness seropositives
  // for(i in 1:500)
  //   seropos_dist_rv[1+500-i] =
  //     1.0 - gamma_cdf(i , seropos_dist_shap, seropos_dist_rate);
  // 
  // // infectiousness
  // pop_infectiousness = conv1d(new_inf, infect_dist_rv);
  // 
  // // seropositives
  // sero_positive = conv1d(new_inf, seropos_dist_rv);
}

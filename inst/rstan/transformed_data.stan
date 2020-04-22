///////////////////////////////////////////////////////////
transformed data {

  int<lower=1>  N_days_tot;
  
  N_days_tot = N_days + N_days_delay; 


  // KEEP COMMENTED OUT: important for later version that take linelist data. 

  // reporting triangle -- laying out cases by date of report and by delay.
  //int           rep_tri_conf_cases[N_days, Max_delay+1]; 
  
  //N_days_tot = N_days + N_days_extra; // how many days to project for
  //for(i in 1:N_days){
  //  for(j in 1:(Max_delay+1)){
  //    rep_tri_conf_cases[i,j] = 0;
  //  }
  //}
  //for(i in 1:N_conf_cases){
  //  rep_tri_conf_cases[cases_test_day[i], cases_days_delay[i]+1] += 1;
  //}
}

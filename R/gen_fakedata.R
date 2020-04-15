genFakeData <- function() { # Should add some params here
  zz <- round(rgamma(1000, 5, 0.3))
  diagnosis_day0 <- max(zz) - zz + 1
  #hist(diagnosis_day0)

  set.seed(124)
  zz <- round(rgamma(1000, 3, 0.9))
  days_delay0 <- zz

  rmv <- (diagnosis_day0 + days_delay0) > max(diagnosis_day0)

  days_delay <- days_delay0[!rmv]
  diagnosis_day <- diagnosis_day0[!rmv]

  # hist(days_delay0); hist(days_delay,add=T,col=5) hist(diagnosis_day0); hist(diagnosis_day,add=T,col=5)
  #range(diagnosis_day)

  # Some real data 
  #raw_data <- read.csv('data/cases_20200325_NYC.csv')
  
   #diagnosis_day0 <- as.numeric(mdy(as.character(raw_data$diagnosis_date)))
   #report_day0    <- as.numeric(mdy(as.character(raw_data$event_create_date)))
   #days_delay     <- report_day0 - diagnosis_day0
   #diagnosis_day  <- diagnosis_day0 - min(diagnosis_day0) + 1 # remove the current day, as seems incomplete
   #rmv            <- max(diagnosis_day) - diagnosis_day - days_delay == 0
   #diagnosis_day  <- diagnosis_day[!rmv]
   #days_delay     <- days_delay[!rmv]

  # Create reporting triangle
  N_days_before <- 10

  # This data is NOT used in the model! It is used to compare to model output
  rep_tri_conf_cases <-
    matrix(0, max(diagnosis_day) + N_days_before, max(days_delay) + 1)

  for (i in 1:length(diagnosis_day)) {
    rep_tri_conf_cases[(diagnosis_day + N_days_before)[i], days_delay[i] + 1] =
      rep_tri_conf_cases[(diagnosis_day + N_days_before)[i], days_delay[i] + 1] + 1
  }

  list(
    diagnosis_day      = diagnosis_day,
    N_days_before      = N_days_before,
    days_delay         = days_delay,
    rep_tri_conf_cases = rep_tri_conf_cases
  )
}


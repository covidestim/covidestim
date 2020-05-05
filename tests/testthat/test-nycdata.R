test_that("NYC data is validated by model", {

  library(tidyr)
  library(dplyr)

  cases <- read.csv("test-nycdata-data.csv", stringsAsFactors = FALSE)

  cases <- replace_na(cases, list("New_COVID_CASE_COUNT" = 0, 
                                  "HOSPITALIZED_CASE_COUNT" = 0,
                                  "DEATH_COUNT"= 0))

  # reformat dates! 
  case <- select(cases, REPORT_DATE, NEW_COVID_CASE_COUNT) %>% 
    rename(date = REPORT_DATE,
    observation = NEW_COVID_CASE_COUNT) %>%
    mutate(date = lubridate::mdy(date)) 

  hosp <- select(cases, REPORT_DATE, HOSPITALIZED_CASE_COUNT) %>% 
          rename(date = REPORT_DATE,
          observation = HOSPITALIZED_CASE_COUNT)%>%
          mutate(date = lubridate::mdy(date)) 
                            
  mort <- select(cases, REPORT_DATE, DEATH_COUNT) %>% 
    rename(date = REPORT_DATE,
    observation = DEATH_COUNT)%>%
    mutate(date = lubridate::mdy(date)) 

  expect_silent(
    cfg <- covidcast(N_days=45, chains=3) +
      input_cases(case) +
      input_hospitalizations(hosp) +
      input_deaths(mort) +
      priors_diagnosis(p_diag_if_inf = c(.1, 9.9), 
                       p_diag_if_sym = c(5,5), 
                       p_diag_if_hos = c(9.5, .5))
  )
})

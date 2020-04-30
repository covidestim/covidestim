#' @export
print.simulateddata <- function(d) {
  cat("True data:\n")
  print(head(d$sim_true))
  cat("\nReported data:\n")
  print(head(d$sim_repor))
}

#' Generate simulated epidemic data 
#'
#' Simulates the first \code{n_days} days of an outbreak of Covid-19, returning
#' reporting data, and the underlying truthy data that the reporting data was
#' generated from.
#'
#' The outbreak starts with 5 new symptomatic cases on day 1 and grows
#' exponentially at a rate of 1.03 – 1.09 to approximately 75,000 total
#' symptomatic cases by day 50, with default parameters. Hospitalizations are
#' simulated as 35% of incident cases, occurring \eqn{x} days after the date of
#' infection, where \eqn{x} is drawn from a gamma distribution gamma (5,0.7),
#' for a total of ~5,000 hospitalizations by day 50. Deaths are simulated as 1%
#' of incident cases, occurring \eqn{y} days after the date of infection, where
#' \eqn{y} is drawn from a gamma distribution \eqn{y \sim \Gamma(\alpha = 4.5,
#' \beta = 0.5)}{gamma (4.5, 0.5)}, for a total of ~1,250 deaths by day 50.
#' ‘reported’ data is created by simulating gaps and delays in detection and
#' reporting of true symptomatic cases. 50\% of symptomatic cases on a given day
#' \eqn{i} will be detected. For each detected case on day \eqn{i}, a date of
#' detection, \eqn{j}, is assigned, and a date of reporting, \eqn{k}, such
#' that:
#'
#' \deqn{\textrm{case}_{\textrm{true},i} = \textrm{case}_{\textrm{reported},
#' i+j+k}}{case_(true,i) = case_(reported, i+j+k)}
#' 
#' Where \eqn{j} and \eqn{k} are drawn from discrete gamma distributions,
#' \eqn{\Gamma(\alpha=4, \beta=1)}{(4,1)} and \eqn{\Gamma(\alpha=2,
#' \beta=1)}{(2,1)} respectively. We simulate hospitalization in a similar
#' manner; 95\% of hospitalizations attributed to a case will be detected. For
#' each detected hospitalization, we assign a date of detection, \eqn{j}, and a
#' date of reporting, \eqn{k}, such that:
#' 
#' \deqn{\textrm{hospital}_{\textrm{true},i} =
#' \textrm{hospital}_{\textrm{reported}, i+j+k}}{case_(true,i) =
#' case_(reported, i+j+k)}
#' 
#' Where \eqn{j} and \eqn{k} are drawn from discrete gamma distributions,
#' \eqn{\Gamma(\alpha=2, \beta=1)}{(2,1)} and \eqn{\Gamma(\alpha=1.7,
#' \beta=0.75)}{(1.7,0.75)} respectively.  Finally, we simulated reported
#' deaths assuming that 95\% of deaths attributed to a will be properly
#' attributed with a reporting delay drawn from a gamma distribution
#' \eqn{\Gamma(\alpha=2, \beta=1)}{(2,1)}. We use simulated data from days 6 –
#' 50 (45 days total) as input into the model.
#'
#' @param p_cases_hospitalized A number in \code{[0,1]}. The fraction of cases
#'   that are hospitalized.
#' @param p_cases_die A number in \code{[0,1]}. The fraction of cases that die.
#' @param p_cases_diagnosed A number in \code{[0,1]}. The fraction of cases that
#'   receive diagnosis outside the hospital.
#' @param p_hospitalizations_diagnosed A number in \code{[0,1]}. The fraction of
#'   hospitalized cases that are diagnosed
#' @param p_deaths_diagnosed A number in \code{[0,1]}. The fraction of deaths that
#'   are diagnosed
#' @param n_days A number. The number of days to simulate data for. The
#'   recommended number is greater than 40.
#'
#' @return A list with two items: \code{reported} and \code{true}.
#'   \itemize{
#'     \item \code{reported} is simulated case reporting data. It is
#'       \code{\link[tibble]{tibble}} containing the variables
#'       \code{reporting_day}, \code{cases}, \code{hosp}, and \code{deaths}.
#'     \item \code{true} is simulated case data, with total knowledge of all
#'       cases. It is also a \code{\link[tibble]{tibble}} containing the
#'       variables \code{reporting_day}, \code{cases}, \code{hosp}, and
#'       \code{deaths}.
#'   }
#'
#'   Note that there are implicit missing values due to some days not being
#'   reported. \code{\link{input_cases}}, \code{\link{input_deaths}}, and 
#'   \code{\link{input_hospitalizations}} all require that missingness in the
#'   data is represented explicitly, using \code{0} to represent a missing
#'   value, and all of these functions will issue an error if there is implicit
#'   missingness. Use a function such as \code{\link[tidyr]{complete}} to fill
#'   in these missing values with zeroes before using output of this function
#'   with Covidcast.
#'
#' @importFrom magrittr %>%
#' @export
simulated_data <- function(p_cases_hospitalized = 0.35,
                           p_cases_die = 0.03,
                           p_cases_diagnosed = 0.5,
                           p_hospitalizations_diagnosed = 0.95,
                           p_deaths_diagnosed = 0.95,
                           n_days = 60) {

  # TRUE symptomatic cases: 
  # number of days of symptomatic cases to simulate
  days <- n_days

  # simulate symptomatic cases (cumulative)
    ccases <- rep(NA, days)
    ccases[1] <- 5

    ## simulate growth rates slowing slightly over time
      for (i in 2:10){
        growth <- runif(10000, 1.07, 1.09) 
        ccases[i] <- ccases[i-1]^(growth[i])
      }
      for(i in 11:20){
        growth <- runif(10000, 1.04, 1.045) 
        ccases[i] <- ccases[i-1]^(growth[i]) 
      }
      for(i in 21:days){
        growth <- runif(10000, 1.025, 1.0275) 
        ccases[i] <- ccases[i-1]^(growth[i]) 
      }
  # convert from cumulative (ccases) to daily count
  cumulative <- as.data.frame(round(ccases)) %>%
                dplyr::mutate(daily = `round(ccases)` - 
                         dplyr::lag(`round(ccases)`, 1)) %>%
                tidyr::replace_na(list(daily = 5))
  # finaly daily case count df, where rowid is day and dialy is count:
  cases_final <- as.data.frame(cumulative$daily) %>% 
                 tibble::rowid_to_column() %>%
                 dplyr::rename(cases = `cumulative$daily`)

  # now we need simulate hospitalizations and deaths
  # based on a probability of progression and a delay

  # first, generate random # and delay daysd

  fate <- runif(cumulative[days,1], 0, 1)
  hosp_delay <- rgamma(length(fate), 5, 0.7)
  die_delay <- rgamma(length(fate), 4.5, 0.5)

  # create a linelist: 
  # each row is a case
  # the value of "ll_case" is the date of symptom onset
  for_ll <- as.data.frame(rep(cases_final$rowid, cases_final$cases)) %>%
            # add random draws
            cbind(fate, hosp_delay, die_delay) %>%
            dplyr::rename(ll_case = `rep(cases_final$rowid, cases_final$cases)`)

  #35% of cases are hospitalized 
  hosp_ll <- dplyr::filter(for_ll, fate <= p_cases_hospitalized) %>% 
             dplyr::mutate(day = round(ll_case + hosp_delay))
  #% of hospitalizations / 2% cases will die
  die_ll <- dplyr::filter(for_ll, fate <= p_cases_die) %>% 
            dplyr::mutate(day = round(ll_case + hosp_delay + die_delay))

  # summarise each _ll back to daily counts: 
  hospital <- dplyr::group_by(hosp_ll, day) %>% dplyr::summarise(hosp = dplyr::n())
  death <- dplyr::group_by(die_ll, day) %>% dplyr::summarise(death = dplyr::n())

  # join together in 'true' dataset
  sim_true <- dplyr::left_join(cases_final, hospital, by = c("rowid" = "day")) %>%
              dplyr::left_join(., death, by = c("rowid" = "day")) %>%
              tidyr::replace_na(list(hosp = 0, death = 0))

  # OBSERVED Data 
  # to be input into package
  # adding detection gaps and delays
  case_dd <- round(rgamma(nrow(for_ll),4,1))
  case_rd <- round(rgamma(nrow(for_ll),2,1))
  hosp_dd <- round(rgamma(nrow(hosp_ll),2,1))
  hosp_rd <- round(rgamma(nrow(hosp_ll),1.7,0.75))
  die_rd  <- round(rgamma(nrow(die_ll),2,1))
  
  

  case_by_report <- cbind(for_ll, case_dd, case_rd) %>%
                     dplyr::filter(fate <= p_cases_diagnosed) %>% #50% of cases are diagnosed
                     dplyr::mutate(report_day = ll_case + case_dd + case_rd) %>%
                     dplyr::group_by(report_day) %>%
                     dplyr::summarise(case = n()) %>%
                     dplyr::filter(report_day <= days)

  hosp_by_report <- cbind(hosp_ll, hosp_dd, hosp_rd) %>%
                    dplyr::filter(fate <= p_hospitalizations_diagnosed) %>% #95% of hospitalizations are diagnosed
                    dplyr::mutate(report_day = day + hosp_dd + hosp_rd) %>%
                    dplyr::group_by(report_day) %>%
                    dplyr::summarise(hospital = n()) %>%
                    dplyr::filter(report_day <= days)

  death_by_report <- cbind(die_ll, die_rd) %>%
                    dplyr::filter(fate <= p_deaths_diagnosed) %>% #95% of deaths are diagnosed
                    dplyr::mutate(report_day = day + die_rd) %>%
                    dplyr::group_by(report_day) %>%
                    dplyr::summarise(death = n()) %>%
                    dplyr::filter(report_day <= days)

  sim_repor <- dplyr::full_join(case_by_report, hosp_by_report, 
                            by = c("report_day" = "report_day")) %>%
                  dplyr::left_join(., death_by_report, 
                            by = c("report_day" = "report_day")) %>%
                  tidyr::replace_na(list(case = 0, 
                                  hospital = 0, 
                                  death = 0)) %>%
                  dplyr::filter(report_day >= 11) %>% 
                  ## NOTE sim_repor starts on day 6 of sim_true
                    # e.g. no reported cases until 6 days after
                    # first truly symptomatic case occurs
                  rbind(c(11, 0 , 0, 0)) %>%
                  dplyr::arrange(report_day) %>%
                  dplyr::rename(cases=case, hosp=hospital, deaths=death)

  sim_true <- dplyr::rename(sim_true, report_day=rowid, cases=cases, deaths=death)
    
  structure(list(reported = sim_repor,
                 true     = tibble::as_tibble(sim_true)),
            class="simulateddata")
}

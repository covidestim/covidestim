#' @export
print.simulateddata <- function(d) {
  cat("True data:\n")
  print(d$true)
  cat("\nReported data:\n")
  print(d$reported)
}

#' Generate simulated epidemic data 
#'
#' Simulates the first \code{n_days} days of an outbreak of Covid-19, returning
#' reporting data, and the underlying truthy data that the reporting data was
#' generated from.
#'
#' The outbreak starts with 10 new symptomatic cases on day 1 and grows
#' at a rate sampled from a log normal distribution such that cases on day i+1
#' are equal to cases on day i multiplied by a growth rate j where j is 
#' distributed log normal (0.175, 0.05). There are a cumulative total of 
#' ~452,000 symptomatic cases by day \code{n_days}, with default parameters. 
#' 
#' Deaths are simulated as 2.5% of incident cases, occurring \eqn{y} days after 
#' the date of symptom onset, where  \eqn{y} is drawn from a gamma distribution 
#' \eqn{y \sim \Gamma(\alpha = 4.4, \beta = 0.22)}{gamma (4.4, 0.22)}, for a 
#' total of ~840 deaths by day \code{n_days}.
#' 
#' ‘Reported’ data is created by simulating gaps and delays in detection and
#' reporting of true symptomatic cases. 50\% of symptomatic cases on a given day
#' \eqn{i} will be detected. For each detected case on day \eqn{i}, a date of
#' detection, \eqn{j}, is assigned, and a date of reporting, \eqn{k}, such
#' that:
#'
#' \deqn{\textrm{case}_{\textrm{true},i} = \textrm{case}_{\textrm{reported},
#' i+j+k}}{case_(true,i) = case_(reported, i+j+k)}
#' 
#' Where \eqn{j} and \eqn{k} are drawn from discrete gamma distributions,
#' \eqn{\Gamma(\alpha=4, \beta=1)}{(4,1)} and \eqn{\Gamma(\alpha=2.4,
#' \beta=0.9)}{(2.4,0.9)} respectively. We simulated reported
#' deaths assuming that 95\% of deaths will be properly
#' attributed with a reporting delay drawn from a gamma distribution
#' \eqn{\Gamma(\alpha=1.9, \beta=1.1)}{(1.9,1.1)}. We use simulated data from days 21 –
#' 65 (45 days total) as input into the model.
#'
#' @param p_cases_die A number in \code{[0,1]}. The fraction of cases that die.
#' @param p_cases_diagnosed A number in \code{[0,1]}. The fraction of cases that
#'   receive diagnosis outside the hospital.
#' @param p_deaths_diagnosed A number in \code{[0,1]}. The fraction of deaths that
#'   are diagnosed
#' @param n_days A number. The number of days to simulate data for. The
#'   recommended number is greater than 40 for simulated case data and over 60
#'   for simulated death data. 
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
simulated_data <- function(p_cases_die = 0.025,
                           p_cases_diagnosed = 0.5,
                           p_deaths_diagnosed = 0.95,
                           n_days = 65) {

  list(
    p_cases_die,
    p_cases_diagnosed,
    p_deaths_diagnosed
  ) -> ps

  # Input validation
  purrr::walk(
    ps,
    ~att(is.numeric(.), . < 1, . > 0)
  )

  # TRUE symptomatic cases: 
  # number of days of symptomatic cases to simulate
  days <- n_days

  # simulate symptomatic cases (cumulative)
    ccases <- rep(NA, days)
    ccases[1] <- 10

    lambda = rlnorm(days - 1, meanlog = 0.175, sdlog = 0.05)
    
    # start loop: Take previous year's N and multiply by lambda
    for (t in 2:days) {
      ccases[t] = ccases[t - 1] * lambda[t - 1]
    }
 
  # convert from cumulative (ccases) to daily count
  cumulative <- as.data.frame(round(ccases)) %>%
                dplyr::mutate(daily = `round(ccases)` - 
                         dplyr::lag(`round(ccases)`, 1)) %>%
                tidyr::replace_na(list(daily = 10))
  # finaly daily case count df, where rowid is day and dialy is count:
  cases_final <- as.data.frame(cumulative$daily) %>% 
                 tibble::rowid_to_column() %>%
                 dplyr::rename(cases = `cumulative$daily`)

  # now we need simulate hospitalizations and deaths
  # based on a probability of progression and a delay

  # first, generate random # and delay daysd

  fate <- runif(cumulative[days,1], 0, 1)
  die_delay <- rgamma(length(fate), 4.4, 0.22)

  # create a linelist: 
  # each row is a case
  # the value of "ll_case" is the date of symptom onset
  for_ll <- as.data.frame(rep(cases_final$rowid, cases_final$cases)) %>%
            # add random draws
            cbind(fate, die_delay) %>%
            dplyr::rename(ll_case = `rep(cases_final$rowid, cases_final$cases)`)

  #% of hospitalizations / 2% cases will die
  die_ll <- dplyr::filter(for_ll, fate <= p_cases_die) %>% 
            dplyr::mutate(day = round(ll_case + die_delay))

  # summarise each _ll back to daily counts: 
  death <- dplyr::group_by(die_ll, day) %>% dplyr::summarise(death = dplyr::n())

  # join together in 'true' dataset
  sim_true <- dplyr::left_join(cases_final, death, by = c("rowid" = "day")) %>%
              tidyr::replace_na(list(hosp = 0, death = 0))

  # OBSERVED Data 
  # to be input into package
  # adding detection gaps and delays
  case_dd <- round(rgamma(nrow(for_ll),4,1))
  case_rd <- round(rgamma(nrow(for_ll),2.4,0.9))
  die_rd  <- round(rgamma(nrow(die_ll),1.9,1.1))
  
  case_by_report <- cbind(for_ll, case_dd, case_rd) %>%
                     dplyr::filter(fate <= p_cases_diagnosed) %>% #50% of cases are diagnosed
                     dplyr::mutate(report_day = ll_case + case_dd + case_rd) %>%
                     dplyr::group_by(report_day) %>%
                     dplyr::summarise(case = n()) %>%
                     dplyr::filter(report_day <= days)

  death_by_report <- cbind(die_ll, die_rd) %>%
                    dplyr::filter(fate <= p_deaths_diagnosed) %>% #95% of deaths are diagnosed
                    dplyr::mutate(report_day = day + die_rd) %>%
                    dplyr::group_by(report_day) %>%
                    dplyr::summarise(death = n()) %>%
                    dplyr::filter(report_day <= days)

  sim_repor <- dplyr::full_join(case_by_report, death_by_report, 
                            by = c("report_day" = "report_day")) %>%
                  tidyr::replace_na(list(case = 0, 
                                  hospital = 0, 
                                  death = 0)) %>%
                  dplyr::filter(report_day >= 21) %>% 
                  ## NOTE sim_repor starts on day 6 of sim_true
                    # e.g. no reported cases until 6 days after
                    # first truly symptomatic case occurs
                  dplyr::arrange(report_day) %>%
                  dplyr::rename(cases=case, deaths=death)

  sim_true <- dplyr::rename(sim_true, report_day=rowid, cases=cases, deaths=death)
    
  structure(list(reported = sim_repor,
                 true     = tibble::as_tibble(sim_true)),
            class="simulateddata")
}

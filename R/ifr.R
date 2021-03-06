#' Estimated IFR ... for US counties and states
#'
#' DETAILS HERE. Look in `data-raw/` for details on how this is generated.
#'
#' @param region A string with the state name, or the FIPS code
#'
#' @return A data.frame with two variables: [date, value], where value is
#'   the odds ratio on VALUE HERE for that particular day. If \code{region}
#'   is specified incorrectly, an error will be thrown
#'
#' @examples
#' get_ifr_raw('Connecticut')
#' get_ifr_raw('09009')
get_ifr_raw <- function(region) {
  found_state <- FALSE
  region_is_county <- any(region %in% ifr_county$fips)
  
  # Branch for case when 'region' is a state
  if (region %in% ifr_state$State)
    found_state <- ifr_state$ifr_OR[ifr_state$State==region]
  # Branch for case when 'region' is not a state and is probably a county
  else if (region_is_county)
    found_state <- local({
      state_name <- dplyr::filter(ifr_county, fips == region)$State[[1]]

      ifr_state$ifr_OR[ifr_state$State==state_name]
    })
  else
    stop("`region` was neither a state name or a character FIPS code")

  successful_state_find <- function(candidate) is.numeric(candidate)

  if (!successful_state_find(found_state))
    stop(glue::glue("Could not find state-level IFR data for region {region}!"))

  if (!region_is_county)
    return(found_state)

  found_county <- dplyr::filter(ifr_county, fips == region)

  if (nrow(found_county) > 1)
    stop(glue::glue("Expected `found_county` to be scalar, but matched more than 1 row"))

  # Check to be sure we found county IFR data
  if (nrow(found_county) == 0)
    stop(glue::glue("Could not find county-level IFR information for county {region}"))

  # Multiply the enclosing state's `value` column by the `comorb_OR` scalar
  found_state * found_county$ifr_OR
}

gen_ifr_adjustments <- function(first_date, N_days_before, region) {
  # The next several lines add state/county-specific IFR data to the
  # configuration passed to Stan. These lines have to be defined inside 
  # this function because they need to know the start date of the input data,
  # which is not known until the user adds `input_cases()/input_deaths()` to
  # the model configuration.
  ymd <- lubridate::ymd

  ifr_adj_start <- first_date - 
    lubridate::days(N_days_before)  # Burn-in period
  
  # time-invariant IFR adjustment accounting for state- and county-level
  # factors, representing an odds-ratio applied to the national average IFR
  ifr_adj_fixed <- get_ifr_raw(region) 

  # reduction in IFR over the course of 2020 due to improvements in care.
  # Operationalized as 1 minus a Normal CDF, with max slope in June 1, 
  # and sd = 45 days. Vector created from ifr_adj_start to end 2022.
  ifr_adj_df <- tibble::tibble(
    date  = seq.Date(ifr_adj_start, ymd('2022-12-31'), by = '1 day'),
    value = 1 - pnorm(
        ifr_adj_start : ymd("2022-12-31"),
        ymd("2020-6-1"),
        45
      )
  )

  ifr_adj <- dplyr::pull(ifr_adj_df, value)

  list(
    ifr_adj_fixed = ifr_adj_fixed,
    ifr_adj       = ifr_adj,
    N_ifr_adj     = length(ifr_adj)
  )
}


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
  if (region %in% colnames(ifr_state))
    found_state <- dplyr::select(ifr_state, date, value = all_of(region))
  # Branch for case when 'region' is not a state and is probably a county
  else if (region_is_county)
    found_state <- local({
      state_name <- dplyr::filter(ifr_county, fips == region)$state[[1]]

      dplyr::select(ifr_state, date, value = all_of(state_name))
    })
  else
    stop("`region` was neither a state name or a character FIPS code")

  successful_state_find <- function(candidate)
    !identical(candidate, FALSE) && ncol(candidate == 2)

  if (!successful_state_find(found_state))
    stop(glue::glue("Could not find state-level IFR data for region {region}!"))

  if (!region_is_county)
    return(found_state)

  found_county <- dplyr::filter(ifr_county, fips == region)

  # Check to be sure we found county IFR data
  if (nrow(found_county) == 0)
    stop(glue::glue("Could not find county-level IFR information for county {region}"))

  # Multiply the enclosing state's `value` column by the `comorb_OR` scalar
  mutate(found_state, value = value * found_county$comorb_OR)
}

get_ifr <- function(region, start_date) {
  all_data <- get_ifr_raw(region)

  min_date <- min(all_data$date)
  max_date <- max(all_data$date)

  att(class(start_date) == 'Date')

  if (!dplyr::between(start_date, min_date, max_date))
    stop(
      glue::glue(
        "start_date was {start_date} but extent of IFR data is {min_date} - {max_date})"
      )
    )

  dplyr::filter(all_data, date >= start_date)
}

gen_ifr_adjustments <- function(first_date, N_days_before, region) {
  # The next several lines add state/county-specific IFR data to the
  # configuration passed to Stan. These lines have to be defined inside 
  # this function because they need to know the start date of the input data,
  # which is not known until the user adds `input_cases()/input_deaths()` to
  # the model configuration.
  ymd <- lubridate::ymd

  ifr_adj_start <- first_date - 
    lubridate::days(N_days_before) + # Burn-in period
    lubridate::days(14)              # 14-days for deaths to start NEEDS BETTER DOCUMENTATION

  # IFR data, beginning at `ifr_adj_start` and ending on 2021-12-31, or the
  # maximum date-value in the CSV files in `data-raw/ifr-data`.
  #
  # The cubed root of these data are taken to represent this adjustment
  # being applied 3 times.
  ifr_adj_df <- get_ifr(region, start_date = ifr_adj_start) %>%
    mutate(value = value^(1/3))

  ifr_adj <- dplyr::pull(ifr_adj_df, value)

  # lower mortality due to improved care (assumes IFR 30% higher
  # at start of pandemic)
  ifr_adj2_df <- tibble::tibble(
    date  = seq.Date(ymd('2020-01-01'), ymd('2021-12-31'), by = '1 day'),
    value = 1 +
      0.30 * (1 - pnorm(
        ymd("2020-1-1") : ymd("2021-12-31"),
        ymd("2020-6-1"),
        45
      ))               
  ) %>% dplyr::inner_join(ifr_adj_df, by='date') %>%
    transmute(date, value = value.x*value.y)

  ifr_adj2 <- dplyr::pull(ifr_adj2_df, value)

  list(
    ifr_adj   = ifr_adj,
    ifr_adj2  = ifr_adj2,
    N_ifr_adj = length(ifr_adj)
  )
}


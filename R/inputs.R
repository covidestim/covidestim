validate_input <- function(d) {
  att(is.data.frame(d))
  att(nrow(d) >= 1)
  att(setequal(names(d), c("date", "observation")))
  att(is.numeric(d$observation))
  att(all(d$observation >= 0))
}

transform_input <- function(d)
  dplyr::mutate(
    d,
    date = reformat_dates(date),
    observation = as.integer(observation)
  )

reformat_dates <- function(vec) lubridate::ymd(vec)

#' Input observational data
#'
#' There are three types of observational data that can be used with Covidcast.
#'
#' \itemize{
#'   \item Case reporting data, detailing the number of new cases each day
#'   \item Hospitalization data, detailing the number of hospitalizations each day
#'   \item Death data, detailing the number of confirmed Covid-19 deaths each day
#' }
#'
#' All input data to Covidcast is expected to be a
#' \code{\link[base]{data.frame}}, of two variables One variable \code{date}
#' will be of the form \code{YYYY-MM-DD}.  The second column,
#' \code{observations} will be a non-negative numeric vector.
#'
#' Missing values should be represented as 0. The date range of the three sets
#' of data must be equivalent, with one observation each day, and no gaps in
#' the data. Several checks will attempt to enforce this specification.
#'
#' @export
input_cases <- function(data) {
  validate_input(data)
  data <- transform_input(data)
  structure(list(obs_cas=data), class='input')
}

#' @rdname input_cases
#' @export
input_deaths <- function(data) {
  validate_input(data)
  data <- transform_input(data)
  structure(list(obs_die=data), class='input')
}

#' @rdname input_cases
#' @export
input_hospitalizations <- function(data) {
  validate_input(data)
  data <- transform_input(data)
  structure(list(obs_hos=data), class='input')
}

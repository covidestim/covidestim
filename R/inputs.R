validate_input <- function(d) {
  att(is.data.frame(d))
  att(nrow(d) >= 1)
  att(setequal(names(d), c("date", "observation")))
  att(is.numeric(d$observation))
  att(all(d$observation >= 0))
}

reformat_dates <- function(vec) lubridate::ymd(vec)

#' Input case data
#'
#' Takes daily diagnosis data as \code{data}.
#'
#' @param data A \code{\link[base]{data.frame}} containing case data
#'
#' @export
input_cases <- function(data) {
  validate_input(data)

  data <- dplyr::mutate(data, date = reformat_dates(date))

  structure(list(obs_cas=data), class='input')
}

#' Input deaths data
#'
#' Takes daily deaths data as \code{data}.
#'
#' @param data A \code{\link[base]{data.frame}} containing deaths data
#'
#' @export
input_deaths <- function(data) {
  validate_input(data)
  
  data <- dplyr::mutate(data, date = reformat_dates(date))

  structure(list(obs_die=data), class='input')
}

#' Takes hospitalizations data as \code{data}.
#'
#' @param data A \code{\link[base]{data.frame}} containing hospitalizations
#'   data
#'
#' @export
input_hospitalizations <- function(data) {
  validate_input(data)
  
  data <- dplyr::mutate(data, date = reformat_dates(date))

  structure(list(obs_hos=data), class='input')
}

#att <- assertthat::assert_that

validate_input <- function(data) {
  att(is.data.frame(data))
  att(nrow(data) > 1)
}

#' Input case data
#'
#' Takes daily diagnosis data as \code{data}.
#'
#' @param data A \code{\link[base]{data.frame}} containing case data
#' @param ddate A string. The variable name containing date data. Date data
#'   should be specified as YYYY-MM-DD, and should be of type \code{character}.
#' @param count A string. The variable name containing the counts of diagnosed
#'   patients. Must be of type \code{numeric}.
#'
#' @export
input_cases <- function(data, ddate='date', count='count') {
  validate_input(data)
  
  data <- dplyr::select(data, ddate=ddate, count=count)

  structure(list(data=data), class='input')
}

#' Input deaths data
#'
#' Takes daily deaths data as \code{data}.
#'
#' @param data A \code{\link[base]{data.frame}} containing deaths data
#' @param ddate A string. The variable name containing date data. Date data
#'   should be specified as YYYY-MM-DD, and should be of type \code{character}.
#' @param count A string. The variable name containing the counts of deaths.
#'   Must be of type \code{numeric}.
#'
#' @export
input_deaths <- function(data, ddate='date', count='count') {
  validate_input(data)
  
  data <- dplyr::select(data, ddate=ddate, count=count)

  structure(list(data=data), class='input')
}

#'
#' Takes hospitalizations data as \code{data}.
#'
#' @param data A \code{\link[base]{data.frame}} containing hospitalizations
#'   data
#' @param ddate A string. The variable name containing date data. Date data
#'   should be specified as YYYY-MM-DD, and should be of type \code{character}.
#' @param count A string. The variable name containing the counts of newly
#'   hospitalized parients. Must be of type \code{numeric}.
#'
#' @export
input_hospitalizations <- function(data, hdate='date', count='count') {
  validate_input(data)
  
  data <- dplyr::select(data, hdate=hdate, count=count)

  structure(list(data=data), class='input')
}

#'
#' Takes cumulative case data as \code{data}.
#'
#' @param data A \code{\link[base]{data.frame}} containing testing data
#' @param ddate A string. The variable name containing date data. Date data
#'   should be specified as YYYY-MM-DD, and should be of type \code{character}.
#' @param count A string. The variable name containing the number of positive
#'   tests recorded that day. Should be of type \code{numeric}
#'
#' @export
input_tests <- function(data, tdate='date', count='count') {
  validate_input(data)
  
  data <- dplyr::select(data, tdate=tdate, count=count)

  structure(list(data=data), class='input')
}


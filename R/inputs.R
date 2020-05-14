validate_input <- function(d, type) {

  pvec <- purrr::partial(paste, ...=, collapse = ', ')

  att(
    is.data.frame(d),
    msg=glue("Input must be a `data.frame`. Your input was a {pvec(class(d))}")
  )
  att(
    nrow(d) >= 1,
    msg="The input data.frame had 0 rows"
  )
  att(
    setequal(names(d), c("date", "observation")),
    msg=glue("The only variables in the data.frame should be 'date' and 'observation'. ",
             "Yours were: `{vars}`", vars = pvec(names(d)))
  )
  att(
    "POSIXct" %in% class(d$date) | "Date" %in% class(d$date),
    msg=glue("The `date` variable must be of class `POSIXct` or `Date`. ",
             "Your `date` variable was of class `{pvec(class(d$date))}`. ",
             "Consider using as.Date()?")
  )
  att(
    is.numeric(d$observation),
    msg=glue(
      "The `observation` variable must be a numeric vector. ",
      "Your `observation` variable was of type ",
      "`{pvec(class(d$observation))}`"
    )
  )
  att(
    all(d$observation >= 0),
    msg=glue("At least one observation was < 0. ",
             "This occurred on rows {pvec(which(d$observation < 0))}")
  )
  att(
    assertthat::is.string(type),
    msg = "Type of input data must be passed as a string."
  )
  att(
    type %in% c("reported", "occurred"),
    msg = glue("`type` of input data must be one of 'reported', 'occurred'. ",
               "You passed {type}.")
  )
}

transform_input <- function(d)
  dplyr::mutate(
    d,
    date        = reformat_dates(date),
    observation = as.integer(observation)
  )

reformat_dates <- function(vec) vec

#' Input observational data
#'
#' There are three types of observational data that can be used with Covidcast.
#'
#' \itemize{
#'   \item Case reporting data, detailing the number of new cases each day
#'   \item Hospitalization data, detailing the number of hospitalizations each
#'   day
#'   \item Death data, detailing the number of confirmed Covid-19 deaths each
#'   day
#' }
#'
#' All input data to Covidcast is expected to be a
#' \code{\link[base]{data.frame}} of two variables. One variable \code{date}
#' must be a vector of type \code{\link[base]{POSIXct}} or
#' \code{\link[base]{Date}}. The second column, \code{observations} will be a
#' non-negative numeric vector.
#'
#' Missing values should be represented as \code{0}. The date range of the three sets
#' of data must be equivalent, with one observation each day, and no gaps in
#' the data. Assertions attempt to enforce this specification.
#'
#' @export
input_cases <- function(data, type = "reported") {
  validate_input(data, type)
  data <- transform_input(data)
  structure(list(obs_cas=data), class='input', date_type = type)
}

#' @rdname input_cases
#' @export
input_deaths <- function(data, type = "reported") {
  validate_input(data, type)
  data <- transform_input(data)
  structure(list(obs_die=data), class='input', date_type = type)
}

#' @rdname input_cases
# input_hospitalizations <- function(data) {
#   validate_input(data)
#   data <- transform_input(data)
#   structure(list(obs_hos=data), class='input')
# }

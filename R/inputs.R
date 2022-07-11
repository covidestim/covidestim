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
#' A family of functions used for inputting data into Covidestim.
#'
#' There are two types of observational data that can be used with Covidestim:
#'
#' \itemize{
#'   \item Case data, detailing the number of new cases each week
#'   \item Death data, detailing the number of confirmed Covid-19 deaths each
#'   week
#'   \item Hospitalizations data, detailing the number of hospital admissions each week
#'   \item Booster data, detailing the number of booster shots adminsitered each week
#'   \item RR data, detailing the Relative Risk IFR adjustment that the model needs to account for.
#' }
#'
#' All input data to Covidestim is expected to be a
#' \code{\link[base]{data.frame}} of two variables. One variable, \code{date},
#' must be a vector of type \code{\link[base]{POSIXct}} or
#' \code{\link[base]{Date}}. The second column, \code{observations}, must be a
#' non-negative numeric vector.
#'
#' Missing values in cases, hospitalizations, booster, RR and deaths data should be imputed before running the model.
#' The date range of all sets of data passed to \code{\link{covidestim}} must be
#' equivalent, with one observation each day, and no gaps in the data.
#' Assertions attempt to enforce this specification.
#'
#' @param data A \code{\link[base]{data.frame}} as described below.
#' @param type A string, either \code{"reported"} or \code{"occurred"}.
#' 
#'   \code{"reported"} applies to the following situations:
#'   \itemize{
#'     \item Case, death, hospitalizatons, booster or RR data where the count for a particular day
#'     represents the number of cases or deaths that were reported publicly on
#'     that day for the following week.
#'     \item Cases, death, hospitalizatons, booster or RR data where the count for a particular day is the 
#'       number of cases or deaths received from hospital reports, or
#'       diagnostic laboratories, accumulated across a particular region
#'   }
#'
#'   \code{"occurred"} applies to the following situations:
#'   \itemize{
#'     \item Data where counts represent the number of cases or deaths that
#'       actually occured on that day. For case data, this would be the day a
#'       test was administered. For deaths, this would be the day an individual
#'       died. Note that, if the deaths data represent the day an individual's
#'       death was first reported as a SARS-Cov-2-related death, that data
#'       should be passed with \code{type = "reported"}, instead.
#'   }
#'
#' @rdname input_cases
#' @export
input_cases <- function(data, type = "reported", lastCaseDate = NULL) {
  validate_input(data, type)
  data <- transform_input(data)
  out <- list(obs_cas=data)

  if(!is.null(lastCaseDate)){
    att(
      "POSIXct" %in% class(lastCaseDate) | "Date" %in% class(lastCaseDate),
      msg=glue("The `lastCaseDate` variable must be NULL or of class `POSIXct` or `Date`. ",
               "Your `lastCaseDate` variable was of class `{class(lastCaseDate)}`. ",
               "Consider using as.Date()?"))
    att(lastCaseDate > data$date[1], msg =glue("The `lastCaseDate` specified was prior to the first observed date. ",
                                                "Specify NULL or a date later than the first observed date"))
    
    attr(out,"lastCaseDate") <- lastCaseDate
  }
  structure(out,class='input', date_type = type)
}

#' @rdname input_cases
#' @export
input_deaths <- function(data, type = "reported",
                         lastDeathDate = NULL) {
  validate_input(data, type)
  data <- transform_input(data)
  out <- list(obs_die=data)
  
  if(!is.null(lastDeathDate)){
    att(
      "POSIXct" %in% class(lastDeathDate) | "Date" %in% class(lastDeathDate),
      msg=glue("The `lastDeathDate` variable must be NULL or of class `POSIXct` or `Date`. ",
               "Your `lastDeathDate` variable was of class `{class(lastDeathDate)}`. ",
               "Consider using as.Date()?"))
    att(lastDeathDate > data$date[1], msg =glue("The `lastDeathDate` specified was prior to the first observed date. ",
                                                "Specify NULL or a date later than the first observed date"))
    
    attr(out,"lastDeathDate") <- lastDeathDate
  }
  structure(out,class='input', date_type = type)
}


#' @rdname input_cases
#' @export
input_rr <- function(data, type = "reported") {
  validate_input(data, type)
  structure(list(ifr_vac_adj=data), class='input', date_type = type)
}


#' @rdname input_cases
#' @export
input_hosp <- function(data, type = "reported",
                       lastHospDate = NULL) {
  
  validate_input(data, type)
  data <- transform_input(data)
  out <- list(obs_hosp=data)
  
  if(!is.null(lastHospDate)){
    att(
      "POSIXct" %in% class(lastHospDate) | "Date" %in% class(lastHospDate),
      msg=glue("The `lastHospDate` variable must be NULL or of class `POSIXct` or `Date`. ",
               "Your `lastHospDate` variable was of class `{class(lastHospDate)}`. ",
               "Consider using as.Date()?"))
    att(lastHospDate > data$date[1], msg =glue("The `lastHospDate` specified was prior to the first observed date. ",
                                                "Specify NULL or a date later than the first observed date"))
    
    attr(out,"lastHospDate") <- lastHospDate
  }
  structure(out,class='input', date_type = type)
}


#' @rdname input_cases
#' @export
input_boost <- function(data, type = "reported") {
  validate_input(data, type)
  structure(list(obs_boost=data), class='input', date_type = type)
}


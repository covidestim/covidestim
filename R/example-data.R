#' Example NYC case/death data
#'
#' Returns a \code{\link[tibble]{tibble}} of NYC case and death
#' data from March 2nd, 2020 to April 15, 2020. Only one type of data will be
#' returned at a time.
#'
#' @param type A string. One of \code{"cases"} or  \code{"deaths"}.
#' 
#' @return A tibble with two variables: \code{date}, a vector of \code{Date} 
#'   objects, and \code{observation}, a vector of doubles.
#'
#' @export
example_nyc_data <- function(type = 'cases') {
  att(type %in% c("cases", "deaths"))

  dplyr::transmute(
    nyc_data, # This is an object from `data-raw/`
    date        = date,
    observation = !!rlang::sym(type) # Hacky stuff (str => symbol)
  )
}

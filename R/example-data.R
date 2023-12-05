#' Example CT case/death data
#'
#' Returns a \code{\link[tibble]{tibble}} of Connecticut case, death, vaccine RR,
#' booster and hospitalizations data from December 2, 2021 to June 30, 2022. 
#' Only one type of data will be returned at a time. 
#'
#' @param type A string: \code{"cases"} or \code{"deaths"} 
#' 
#' @return A tibble with two variables: \code{date}, a vector of \code{Date} 
#'   objects, and \code{observation}, a vector of doubles.
#'
#' @examples
#'
#' cfg <- covidestim(nweeks = 31, region = 'Connecticut',
#'    pop = get_pop("Connecticut")) +
#'   input_cases(example_ct_data('cases')) +
#'   input_deaths(example_ct_data('deaths')) 
#'
#' \dontrun{
#'   result <- run(cfg)
#' }
#'
#' @export
example_ct_data <- function(type = 'cases') {
  att(type %in% c("cases", "deaths"))

  dplyr::transmute(
    ct_data, # This is an object from `data-raw/`
    date        = date,
    observation = !!rlang::sym(type) # Hacky stuff (str => symbol)
  )
}

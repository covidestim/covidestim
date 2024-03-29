#' Example CT case/death data
#'
#' Returns a \code{\link[tibble]{tibble}} of Connecticut case, death, vaccine RR,
#' booster and hospitalizations data from December 2, 2021 to June 30, 2022. 
#' Only one type of data will be returned at a time. 
#'
#' @param type A string: \code{"cases"} or \code{"deaths"} or \code{"RR"} or
#'    \code{"hosp"} or \code{"boost"}
#' 
#' @return A tibble with two variables: \code{date}, a vector of \code{Date} 
#'   objects, and \code{observation}, a vector of doubles.
#'
#' @examples
#'
#' cfg <- covidestim(nweeks = 31, region = 'Connecticut',
#'    pop = get_pop("Connecticut"),
#'    start_p_imm = get_imm_init("Connecticut")$start_p_imm,
#'    cum_p_inf_init = get_imm_init("Connecticut")$cum_p_inf_init) +
#'   input_cases(example_ct_data('cases')) +
#'   input_deaths(example_ct_data('deaths')) +
#'   input_rr(example_ct_data('RR')) + 
#'   input_hosp(example_ct_data('hosp')) +
#'   input_boost(example_ct_data('boost'))
#'
#' \dontrun{
#'   result <- run(cfg)
#' }
#'
#' @export
example_ct_data <- function(type = 'cases') {
  att(type %in% c("cases", "deaths", "RR", "boost", "hosp"))

  dplyr::transmute(
    ct_data, # This is an object from `data-raw/`
    date        = date,
    observation = !!rlang::sym(type) # Hacky stuff (str => symbol)
  )
}

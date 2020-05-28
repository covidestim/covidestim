#' @export
summary.covidcast_result <- function(ccr, include.before = FALSE, 
                                     include.RtEstim = TRUE) {

  att(identical(include.before, FALSE))
  att(identical(include.RtEstim, FALSE))

  c(
    "incidence.est" = "new_inf",
    "cases.fitted"  = "occur_cas",
    "deaths.fitted" = "occur_die",
  ) -> params

  rstan::summary
}

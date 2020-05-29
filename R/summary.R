#' @export
summary.covidcast_result <- function(ccr, include.before = FALSE, 
                                     include.RtEstim = TRUE) {

  att(identical(include.before, FALSE))
  att(identical(include.RtEstim, TRUE))

  c(
    "incidence.est" = "new_inf",
    "cases.fitted"  = "occur_cas",
    "deaths.fitted" = "occur_die"
  ) -> params

  rstan::summary(
    ccr$result,
    pars = "new_inf[1]",
    probs = c(0.025, 0.5, 0.975)
  )
}

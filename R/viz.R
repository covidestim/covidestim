viz <- function(...) UseMethod('viz')

viz.covidcast_result <- function(cc) {

  pdf()

  viz_comparison_to_data_2(cc)
  viz_comparison_to_data_2(cc)
  viz_current_and_cumulative_outcomes(cc)
  viz_incident_outcomes(cc)
  viz_priors_vs_posteriors(cc)

  dev.off()
}

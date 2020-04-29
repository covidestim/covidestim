viz <- function(...) UseMethod('viz')

viz.covidcast_result <- function(cc) {

  pdf()

  viz_melanie(cc)

  dev.off()
}

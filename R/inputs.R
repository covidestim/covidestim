#att <- assertthat::assert_that

validate_input <- function(data) {
  att(is.data.frame(data))
  att(nrow(data) > 1)
}

input_cases <- function(data, ddate='date', count='count') {
  validate_input(data)
  
  data <- dplyr::select(data, ddate=ddate, count=count)

  structure(list(data=data), class='input')
}

input_deaths <- function(data, ddate='date', count='count') {
  validate_input(data)
  
  data <- dplyr::select(data, ddate=ddate, count=count)

  structure(list(data=data), class='input')
}

input_hospitalizations <- function(data, hdate='date', count='count') {
  validate_input(data)
  
  data <- dplyr::select(data, hdate=hdate, count=count)

  structure(list(data=data), class='input')
}

input_tests <- function(data, tdate='date', count='count') {
  validate_input(data)
  
  data <- dplyr::select(data, tdate=tdate, count=count)

  structure(list(data=data), class='input')
}


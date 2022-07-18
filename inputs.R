validate_df <- function(df){
  
  pvec <- purrr::partial(paste, ...=, collapse = ',')
  
  #General data frame test cases 
  att(
    is.data.frame(df),
    msg=glue("Input must be a `data.frame`. Your input was a {pvec(class(df))}")
  )
  
  att(
    nrow(df) >= 1,
    msg="The input data.frame had 0 rows"
  )
  
  att(
    all(df >= 0),
    msg=glue("At least one observation was < 0. ",
             "This occurred on rows {pvec(which(df < 0))}")
  )
  
  att(
    setequal(names(df), c("date", "case", "hospi", "deaths", "boost", "RR")),
    msg=glue("The only variables in the data.frame should be 'date', 'case', 'hospi', 'deaths', 'boost', 'RR'. ",
             "Yours were: `{vars}`", vars = pvec(names(d))
    )
  )
  
  
  #Convert data frame into individual df's to test 
  #Convert date to date type
  date = as.Date(df$date)
  
  #Convert indivdual df's to numeric vectors
  obs_cas = (df[['case']])
  obs_hosp = (df[['hospi']])
  obs_die = (df[['deaths']])
  obs_boost = (df[['boost']])
  ifr_vac_adj = (df[['RR']])
  
  #Test cases for Date 
  att(
    "POSIXct" %in% class(date) | "Date" %in% class(date),
    msg=glue("The `date` variable must be of class `POSIXct` or `Date`. ",
             "Your `date` variable was of class `{pvec(class(date))}`. ",
             "Consider using as.Date()?")
  )
  
  #Test cases for Case 
  att(
    is.numeric(obs_cas),
    msg=glue(
      "The `observation` variable must be a numeric vector. ",
      "Your `observation` variable was of type ",
      "`{pvec(class(obs_cas))}`"
    )
  )
  
  #Test cases for Hospi 
  att(
    is.numeric(obs_hosp),
    msg=glue(
      "The `observation` variable must be a numeric vector. ",
      "Your `observation` variable was of type ",
      "`{pvec(class(obs_hosp))}`"
    )
  )
  
  #Test cases for Deaths 
  att(
    is.numeric(obs_die),
    msg=glue(
      "The `observation` variable must be a numeric vector. ",
      "Your `observation` variable was of type ",
      "`{pvec(class(obs_die))}`"
    )
  )
  
  #Test cases for Boost
  att(
    is.numeric(obs_boost),
    msg=glue(
      "The `observation` variable must be a numeric vector. ",
      "Your `observation` variable was of type ",
      "`{pvec(class(obs_boost))}`"
    )
  )
  
  #Test cases for RR
  att(
    is.numeric(ifr_vac_adj),
    msg=glue(
      "The `RR` variable must be a numeric vector. ",
      "Your `RR` variable was of type ",
      "`{pvec(class(ifr_vac_adj))}`"
    )
  )
}

# validate_input <- function(d, type) {
# 
#   pvec <- purrr::partial(paste, ...=, collapse = ', ')
# 
#   att(
#     is.data.frame(d),
#     msg=glue("Input must be a `data.frame`. Your input was a {pvec(class(d))}")
#   )
#   att(
#     nrow(d) >= 1,
#     msg="The input data.frame had 0 rows"
#   )
#   att(
#     setequal(names(d), c("date", "observation")),
#     msg=glue("The only variables in the data.frame should be 'date' and 'observation'. ",
#              "Yours were: `{vars}`", vars = pvec(names(d)))
#   )
#   att(
#     "POSIXct" %in% class(d$date) | "Date" %in% class(d$date),
#     msg=glue("The `date` variable must be of class `POSIXct` or `Date`. ",
#              "Your `date` variable was of class `{pvec(class(d$date))}`. ",
#              "Consider using as.Date()?")
#   )
#   att(
#     is.numeric(d$observation),
#     msg=glue(
#       "The `observation` variable must be a numeric vector. ",
#       "Your `observation` variable was of type ",
#       "`{pvec(class(d$observation))}`"
#     )
#   )
#   att(
#     all(d$observation >= 0),
#     msg=glue("At least one observation was < 0. ",
#              "This occurred on rows {pvec(which(d$observation < 0))}")
#   )
#   att(
#     assertthat::is.string(type),
#     msg = "Type of input data must be passed as a string."
#   )
#   att(
#     type %in% c("reported", "occurred"),
#     msg = glue("`type` of input data must be one of 'reported', 'occurred'. ",
#                "You passed {type}.")
#   )
# }

transform_input <- function(d)
  dplyr::mutate(
    d,
    date        = reformat_dates(date),
    observation = as.integer(observation)
  )

reformat_dates <- function(vec) vec



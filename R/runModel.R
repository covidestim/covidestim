#' High level description of the function
#'
#' More extended description of the function
#'
#' @param param1 Description
#'
#' @param param2 Description
#'
#' @param param3 Description
#'
#' @return The return value
#'
#' @examples
#' print(mtcars)
#' @export
runModel <- function(data = defaultData(), chains=3, iter=500) {
  ###### ###### SETTINGS ###### ###### ######
  rstan::rstan_options(auto_write = TRUE) # save the compiled executable to '.'
  options(mc.cores = parallel::detectCores())

  ###### ###### RUN MODEL ###### ###### ######

#fit_stan <-
  rstan::stan(
    file    = system.file("rstan/covid_stan_script_MHC_V2.stan",
                         package="covidcast",
                         mustWork=TRUE),
    control = list(adapt_delta = 0.92, max_treedepth = 12),
    data    = data,
    seed    = 1234,
    chains  = chains,
    iter    = iter,
    warmup  = round(0.8*iter)
  )
}

# Example model configuration object. $c can be modified through overloading
# the addition operator
modelconf <- function(a=1, b=2, c=list()) {
  initialArgs <- list(a=a,b=b,c=c)

  # Declare that this is a class
  structure(initialArgs, class='modelconf')
}

# A list with arbitrary values, of class 'argpkg'
argpkg <- function(...) {
  structure(list(...), class='argpkg')
}

# An overloaded addition operator, dispatched on the type of the lhs argument.
# Flips the order of the arguments and calls 'demo_add' to enable type matching
# on the rhs operand.
"+.modelconf" <- function(a, b) {
  # Dispatching on the type of 'b', which should be 'argpkg' in this example
  demo_add(b, a)
}

demo_add <- function(rightside, leftside) UseMethod('demo_add') 

# Specialization for 'argpkg' classes. Splices in all the keys from 'argpkg'
# obj into existing $c of the 'modelconf' object.
demo_add.argpkg <- function(rightside, leftside) {
  leftside$c <- rlang::dots_list(!!!leftside$c, !!!rightside, .homonyms='last')

  leftside
}          

# samps <- rstan::extract(fit_stan)
# 
# dev.off()
# 
# system(paste("open", pdfnam))

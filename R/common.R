# Alias for making assertions on code
att <- assertthat::assert_that

att_w <- function(assertion, msg) {
  tryCatch(
    att(assertion),
    error = function(c) warning(msg)
  )
}

# Combines 'defaults' and 'args' together. If 'args'==list(), it is stripped
# from the result. Any duplicate keys are resolved by selecting the last key
# specified. The result has class 'modelconfig'
splice_class <- function(defaults, args, class)
  structure(
    rlang::dots_list(
      !!! defaults,
      !!! args,
      .homonyms     = 'last',
      .ignore_empty = 'trailing'),
    class=class)


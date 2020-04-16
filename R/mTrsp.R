# Convert a color into an rgb value? TODO: what does this do?
mTrsp <- function(cl, a) {
  apply(grDevices::col2rgb(cl), 2, function(x) {
    rgb(x[1], x[2], x[3], a, maxColorValue = 255)
  })
}


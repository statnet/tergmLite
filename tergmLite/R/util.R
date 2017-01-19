
deinf <- function (x, replace = 1/.Machine$double.eps) {
  ifelse(is.nan(x) | abs(x) < replace, x, sign(x) * replace)
}

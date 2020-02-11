#' code posted by Marat Talipov in stackoverflow post link: https://stackoverflow.com/questions/34445106/test-if-vector-is-contained-in-another-vector-including-repetitions
#'@keywords internal
"%contain%" <- function(values,x) {
  tx <- table(x)
  tv <- table(values)
  z <- tv[names(tx)] - tx
  all(z >= 0 & !is.na(z))
}
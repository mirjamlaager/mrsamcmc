#' @export
logit <- function(p) {
  return (log(p / (1 - p)))
}

#' @export
invlogit <- function(p) {
  return (1 / (1 + exp(-p)))
}

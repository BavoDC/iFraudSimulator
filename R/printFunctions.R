#' Print method for object of type BiRankFr
#'
#' @param x Object of type BiRankFr
#' @param ... Ellipses argument, to be passed to \code{print} function.
#'
#' @method print BiRankFr
#' @return Nothing
print.BiRankFr <- function(x, ...) {
  cat("\n\nResults claims\n\n")
  print(head(x$ResultsClaims), ...)
  cat("\n\nResults parties\n\n")
  print(head(x$ResultsParties), ...)
}

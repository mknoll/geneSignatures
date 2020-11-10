#' @title Outsect
#' @description Outersect
#' with two columns.#' This function does the opposite of intersect()
#' @export
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

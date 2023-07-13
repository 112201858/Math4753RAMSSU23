#' @title My first function
#'
#' @param x A quantitative vector
#'
#' @return A quantitative vector, the square of the argument passed as parameter
#' @export
#'
#' @examples myfun(x = 1:10)
myfun <- function(x) {
  x^2
}




# We don't always want to export a function
  # Sometimes functions are helper functions for the package developer
    # These don't need to be exported

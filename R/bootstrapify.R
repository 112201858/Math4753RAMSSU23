#' Bootstrapify
#'
#' @param iter number of iterations
#' @param x vector from which the function samples
#' @param fun function used for point estimate
#' @param alpha value for confidence interval
#' @param cx scaling factor for graphics
#' @param ... additional parameters
#'
#' @return an invisible list containing {x, a confidence interval, and the function} and a histogram
#' @export
#'
#' @examples
#' bootstrapify(iter=10000,x = 1:100,fun="mean",alpha=0.05,cx=1.5)
bootstrapify <- function(iter = 10000, x, fun = "mean", alpha = 0.05, cx = 1.5, ...){

  n = length(x)

  y = sample(x, n*iter, replace = TRUE)

  rs.mat = matrix(y, nrow = n, ncol = iter, byrow = TRUE)

  xstat = apply(rs.mat, 2, fun)

  ci = quantile(xstat, c(alpha / 2, 1 - alpha / 2))

  para = hist(xstat,
              freq = FALSE,
              las = 1,
              col = "cyan",
              main = paste("Histogram of Bootstrap sample statistics","\n","alpha=",alpha,"iter=",iter,sep=""),
              ...)


  mat = matrix(x, nrow = length(x), ncol = 1, byrow = TRUE)

  pte = apply(mat, 2, fun)

  abline(v = pte, lwd = 3, col = "Black")

  segments(ci[1], 0, ci[2], 0, lwd = 4)

  text(ci[1], 0.025, paste("(", round(ci[1], 2), sep=""), col = "red", cex = cx)

  text(ci[2], 0.025, paste(round(ci[2], 2) , ")", sep=""), col = "red",cex = cx)

  text(pte, max(para$density) / 2, round(pte, 2), cex = cx)

  invisible(list(ci = ci, fun = fun, x = x))

}

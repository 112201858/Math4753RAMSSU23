#' Piecewise Linear
#'
#' @param x vector of data
#' @param coef vector of beta hats
#' @param k knot
#'
#' @return nothing
#' @export
#'
#' @examples
#' piecewiseLinear(1:10,c(1,2,3),4)
piecewiseLinear = function(x,coef,k){

  coef[1]+coef[2]*(x) + coef[3]*(x-k)*(x-k>0)

}

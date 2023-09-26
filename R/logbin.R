#' logbin
#'
#' @param x successes
#' @param param parameter of interest (p in this case)
#' @param n number of trials
#'
#' @return none
#' @export
#'
#' @examples
#' logbin()
logbin=function(x = 5, param = seq(0,1,length=1000), n = 20){

  log(dbinom(x,prob=param,size=n))

}

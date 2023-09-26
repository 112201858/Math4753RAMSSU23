#' My Normal Curve
#'
#' @param mu - mean of a normal distribution
#' @param sigma - standard deviation of a normal distribution
#' @param a - P(Y <= a)
#'
#' @return curve with area shaded and a probability
#' @export
#'
#' @examples
#' myncurve()
myncurve = function(mu = 0, sigma = 1, a = 0){

  curve(dnorm(x,mean=mu,sd=sigma),
        xlim = c(mu-3*sigma, mu + 3*sigma))

  xcurve = seq(mu-3*sigma,
               a,
               length=1000)

  ycurve = dnorm(xcurve,
                 mean=0,
                 sd=1)

  polygon(c(mu-3*sigma,xcurve,a),
          c(0,ycurve,0),col="Red")

  prob = pnorm(a,mu,sigma)
  prob
}

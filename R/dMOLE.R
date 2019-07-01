#' The Marshall-Olkin logistic-exponential distribution
#'
#' @description
#' Density, distribution function, quantile function,
#' random generation and hazard function for the Marshall-Olkin logistic-exponential distribution
#' with parameters \code{mu}, \code{sigma}, and \code{nu} .
#'
#' @param x,q	vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations.
#' @param mu parameter.
#' @param sigma parameter.
#' @param nu parameter.
#' @param log,log.p	logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#'
#' @details
#' The Marshall-Olkin logistic-exponential Distribution with parameters \code{mu},
#' \code{sigma}, and \code{nu} has density given by
#'
#' \eqn{f(x)=\ \mu \sigma \nu \exp^{\nu x}[\exp^{\nu x}-1]^{-\mu-1} / [1 + \sigma [\exp^{\nu x}-1]^{-\mu}]^2}
#'
#' for x > 0
#'
#'
#' @return
#' \code{dMOLE} gives the density, \code{pMOLE} gives the distribution
#' function, \code{qMOLE} gives the quantile function, \code{rMOLE}
#' generates random deviates and \code{hMOLE} gives the hazard function.
#'
#' @examples
#'
#' ## The probability density function
#' par(mfrow=c(1, 2))
#'  curve(dMOLE(x, mu=0.6, sigma=5.5, nu=1), from=0, to=8,
#'  ylim=c(0, 0.3), col="red", las=1, ylab="f(x)")
#'
#'  curve(dMOLE(x, mu=3, sigma=15, nu=1.2), from=0, to=3,
#'  ylim=c(0, 1.5), col="red", las=1, ylab="f(x)")
#'
#' ## The cumulative distribution and the Reliability function
#' par(mfrow=c(1, 2))
#' curve(pMOLE(x, mu=0.6, sigma=5.5, nu=1),
#'       from=0, to=15,  col="red", las=1, ylab="F(x)")
#' curve(pMOLE(x, mu=0.6, sigma=5.5, nu=1, lower.tail=FALSE),
#'       from=0, to=15, col="red", las=1, ylab="S(x)")
#'
#' ## The quantile function
#' p <- seq(from=0, to=0.99999, length.out=100)
#' plot(x=qMOLE(p, mu=0.6, sigma=5.5, nu=1), y=p, xlab="Quantile",
#'      las=1, ylab="Probability")
#' curve(pMOLE(x, mu=0.6, sigma=5.5, nu=1), from=0, add=TRUE, col="red")
#'
#' ## The random function
#' hist(rMOLE(n=10000, mu=0.6, sigma=5.5, nu=1), freq=FALSE,
#'      xlab="x", las=1, main="")
#' curve(dMOLE(x,  mu=0.6, sigma=5.5, nu=1), from=0, add=TRUE, col="red")
#'
#' ## The Hazard function
#' par(mfrow=c(1,1))
#' curve(hMOLE(x,  mu=0.6, sigma=5.5, nu=1), from=0, to=8,
#'       col="red", ylab="Hazard function", las=1)
#'
#' @export
dMOLE <- function(x, mu, sigma, nu, log=FALSE){
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))

  loglik <- log (mu*sigma*nu*exp(nu*x)*(exp(nu*x) -1)^(-mu-1)) - log((1+sigma*(exp(nu*x)-1)^-mu)^2)
  if (log == FALSE)
    density <- exp(loglik)
  else
    density <- loglik
  return(density)
}
#' @export
#' @rdname dMOLE
pMOLE <- function(q, mu, sigma, nu,
                  lower.tail=TRUE, log.p=FALSE){
  if (any(q < 0))
    stop(paste("q must be positive", "\n", ""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))

  cdf <-(1 + sigma*((exp(nu*q)-1)^-mu))^-1

  if (lower.tail == TRUE) cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
#' @export
#' @rdname dMOLE
qMOLE <- function(p, mu, sigma, nu,
                  lower.tail=TRUE, log.p=FALSE){
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))

  if (log.p == TRUE) p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))

  q <- nu^-1 * (log((1 + ((sigma*p)/(1-p)))^(1/mu)))
  q
}
#' @importFrom stats runif
#' @export
#' @rdname dMOLE
rMOLE <- function(n, mu, sigma, nu){
  if(any(n <= 0))
    stop(paste("n must be positive","\n",""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))

  n <- ceiling(n)
  p <- runif(n)
  r <- qMOLE(p, mu, sigma, nu)
  r
}
#' @export
#' @rdname dMOLE
hMOLE <- function(x, mu, sigma, nu){
  if (any(x < 0))
    stop(paste("x must be positive", "\n", ""))
  if (any(mu <= 0 ))
    stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))
    stop(paste("sigma must be positive", "\n", ""))
  if (any(nu <= 0))
    stop(paste("nu must be positive", "\n", ""))

  h <- mu*nu*exp(nu*x) / ((exp(nu*x)-1)*(1 + sigma*(exp(nu*x)-1)^-mu))
  h
}

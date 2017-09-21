\name{SSLASSO}
\alias{SSLASSO}
\title{The Spike-and-Slab LASSO}
\description{Spike-and-Slab LASSO is a spike-and-slab refinement of the LASSO procedure, using a mixture of Laplace priors indexed by  \code{lambda0} (spike) and \code{lambda1} (slab). 
	
	The \code{SSLASSO} procedure fits coefficients paths for Spike-and-Slab LASSO-penalized
  linear regression models over a grid of values for the regularization
  parameter \code{lambda_0}. The code has been adapted from the \code{ncvreg} package (Breheny and Huang, 2011). }
\usage{
SSLASSO(X, y, penalty = c("separable", "adaptive"), variance = c("known", "unknown"), 
lambda1, lambda0, theta = 0.5, sigma = 1, a = 1, b = p, nu = 0, xi = 0, 
eps = 0.001, max.iter = 500,  counter = 10, warn = FALSE)
}
\arguments{
  \item{X}{The design matrix (n x p), without an intercept.  \code{SSLASSO}
  standardizes the data  by default (mean zero and variance \code{n}).}
  \item{y}{Vector of continuous responses (n x 1). The  responses will be centered by default.}
  \item{penalty}{The penalty to be applied to the model.  Either "separable"
    (with a fixed \code{theta}),  "adaptive" (with a random \code{theta}, where \code{theta~ B(a,p)}).}
  \item{variance}{Whether the error variance is also estimated. Either "known" (with a fixed \code{sigma}) or "unknown" (with a random \code{sigma}, where \code{sigma ~ IG(nu, xi)}).}
  \item{lambda1}{Slab variance parameter. Needs to be greater than \code{lambda0}.}
  \item{lambda0}{Spike penalty parameters (L x 1). Either a numeric value for a single run (L=1) or a sequence of increasing values for dynamic posterior exploration.}
  \item{theta}{Prior mixing proportion. For "separable" penalty, this value is fixed. For "adaptive" penalty, this value is used as a starting value. }
  \item{sigma}{Error variance. For "known" variance, this value is fixed. For "unknown" variance, this value is used as a starting value.}
  \item{a}{Hyperparameter of the beta prior \code{B(a,b)} for the adaptive penalty (default \code{a=1}).}
  \item{b}{Hyperparameter of the beta prior \code{B(a,b)} for the adaptive penalty (default \code{b=p}).}
  \item{nu}{Hyperparameter of the inverse-gamma prior \code{IG(nu, xi)} for the unknown variance case (default \code{nu = 0}).}
  \item{xi}{Hyperparameter of the inverse-gamma prior \code{IG(nu, xi)} for the unknown variance case (default \code{xi = 0}).}
  \item{eps}{Optional. Convergence criterion: converged when difference in regression coefficients is less than \code{eps} (default \code{eps = 0.001}).}
  \item{max.iter}{Optional. Maximum number of iterations.  Default is 1000.}
  \item{counter}{Optional. Applicable only for the adaptive penalty. Determines how often the parameter \code{theta} is updated throughout the cycles of coordinate ascent. Default is 10.}
  \item{warn}{TRUE if warnings should be printed; FALSE by default}
}
\details{
  The sequence of models indexed by the regularization parameter
  \code{lambda_0} is fitted using a coordinate descent algorithm. The algorithm uses 
  screening rules for discarding irrelevant predictors along the lines of Breheny (2011).
  
   }
\value{
  An object with S3 class \code{"SSLASSO"} containing:
  \item{beta}{The fitted matrix of coefficients (p x L).  The number of rows is
    equal to the number of coefficients \code{p}, and the number of columns is
    equal to \code{L} (the length of \code{lambda_0}).}
  \item{iter}{A vector of length \code{L} containing the number
    of iterations until convergence at each value of \code{lambda_0}.}
  \item{lambda0}{The sequence of regularization parameter values in the
    path.}
  \item{penalty}{Same as above.}
  \item{thetas}{A vector of length \code{L} containing the hyper-parameter values \code{theta} (the same as \code{theta} for "separable" penalty)}
  \item{sigmas}{A vector of length \code{L} containing the values \code{sigma} (the same as \code{sigma} for "known" variance)}
  \item{select}{A (p x L) binary matrix indicating which variables were selected along the solution path}
  \item{model}{A single model chosen after the stabilization of the regularization path.}
 }
\references{Rockova, V. and George, E.I. (2016) The Spike-and-Slab LASSO. Journal of the American Statistical Association (Accepted).

Moran, G., Rockova, V. and George, E.I. (2017) On variance estimation for high dimensional Bayesian regression. Manuscript.}
\author{Veronika Rockova <Veronika.Rockova@chicagobooth.edu>}
\seealso{\code{\link{plot.SSLASSO}}}
\examples{
## Linear regression, where p > n

library(SSLASSO)

p <- 1000
n <- 100

lambda1 <- 1     # slab penalty for Spike-and-Slab LASSO
lambda0 <- seq(lambda1, 50, length=20) # slab penalties for Spike-and-Slab LASSO
L <- length(lambda0)
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
beta <- c(1, 2, 3, rep(0, p-3))
y = X[,1] * beta[1] + X[,2] * beta[2] + X[,3] * beta[3] + rnorm(n)

# Oracle SSLASSO with known variance
result1 <- SSLASSO(X, y, penalty = "separable", variance = "known", 
               lambda1 = lambda1, lambda0 = lambda0, 
               theta = 3/p)
plot(result1)

# Adaptive SSLASSO with known variance

result2 <- SSLASSO(X, y, penalty = "adaptive", 
               lambda1 = lambda1, lambda0 = lambda0, 
               theta = 0.5, a = 1, b = p, counter = 10)
plot(result2)

# Adaptive SSLASSO with unknown variance
result3 <- SSLASSO(X, y, penalty = "adaptive", variance = "unknown", 
               lambda1 = lambda1, lambda0 = lambda0)
plot(result3)
result3$sigmas[L]
}

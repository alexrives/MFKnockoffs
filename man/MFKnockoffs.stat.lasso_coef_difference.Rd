% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stats_lasso_cv.R
\name{MFKnockoffs.stat.lasso_coef_difference}
\alias{MFKnockoffs.stat.lasso_coef_difference}
\title{Cross-validated penalized linear regression statistics for MFKnockoffs}
\usage{
MFKnockoffs.stat.lasso_coef_difference(X, X_k, y, cores = 2, ...)
}
\arguments{
\item{X}{original design matrix (size n-by-p)}

\item{X_k}{knockoff matrix (size n-by-p)}

\item{y}{response vector (length n). It should be numeric}

\item{cores}{Number of cores used to compute the knockoff statistics by running cv.glmnet.
If not specified, the number of cores is set to approximately half of the number of cores 
detected by the parallel package.}

\item{...}{additional arguments specific to 'glmnet' (see Details)}
}
\value{
A vector of statistics \eqn{W} (length p)
}
\description{
Fit a linear regression model via penalized maximum likelihood and cross-validation.
Then, compute the difference statistic
  \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the 
jth variable and its knockoff, respectively. The value of the regularization
parameter \eqn{\lambda} is selected by cross-validation and computed with glmnet.
}
\details{
This function uses the \code{glmnet} package to fit the lasso path.

This function is a wrapper around the more general \link{MFKnockoffs.stat.glmnet_coef_difference}.

The knockoff statistics \eqn{W_j} are constructed by taking the difference 
between the coefficient of the j-th variable and its knockoff.
 
By default, the value of the regularization parameter is chosen by 10-fold cross-validation.

The optional \code{nlambda} parameter can be used to control the granularity of the 
grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{100},
where \code{p} is the number of columns of \code{X}.

Unless a lambda sequence is provided by the user, this function generates it on a 
log-linear scale before calling 'glmnet' (default 'nlambda': 100).

For a complete list of the available additional arguments, see \link[glmnet]{cv.glmnet}
and \link[glmnet]{glmnet}.
}
\examples{
p=100; n=200; k=15
mu = rep(0,p); Sigma = diag(p)
X = matrix(rnorm(n*p),n)
nonzero = sample(p, k)
beta = 3.5 * (1:p \%in\% nonzero)
y = X \%*\% beta + rnorm(n)

knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
# Basic usage with default arguments
result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, 
                           statistic=MFKnockoffs.stat.lasso_coef_difference)
print(result$selected)

# Advanced usage with custom arguments
foo = MFKnockoffs.stat.lasso_coef_difference
k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
print(result$selected)

}
\seealso{
Other statistics for knockoffs: \code{\link{MFKnockoffs.stat.forward_selection}},
  \code{\link{MFKnockoffs.stat.glmnet_coef_difference}},
  \code{\link{MFKnockoffs.stat.glmnet_lambda_difference}},
  \code{\link{MFKnockoffs.stat.lasso_coef_difference_bin}},
  \code{\link{MFKnockoffs.stat.lasso_lambda_difference_bin}},
  \code{\link{MFKnockoffs.stat.lasso_lambda_difference}},
  \code{\link{MFKnockoffs.stat.random_forest}},
  \code{\link{MFKnockoffs.stat.sqrt_lasso}},
  \code{\link{MFKnockoffs.stat.stability_selection}}
}

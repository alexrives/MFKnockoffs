#' Cross-validated GLM statistics for MFKnockoffs
#' 
#' Fit a generalized linear model via penalized maximum likelihood and cross-validation.
#' Then, compute the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the coefficient estimates for the 
#' jth variable and its knockoff, respectively. The value of the regularization
#' parameter \eqn{\lambda} is selected by cross-validation and computed with glmnet.
#' 
#' @param X original design matrix (size n-by-p)
#' @param X_k knockoff matrix (size n-by-p)
#' @param y response vector (length n). Quantitative for family="gaussian", 
#' or family="poisson" (non-negative counts). For family="binomial" 
#' should be either a factor with two levels, or a two-column matrix of counts 
#' or proportions (the second column is treated as the target class; for a factor, 
#' the last level in alphabetical order is the target class). For family="multinomial", 
#' can be a nc>=2 level factor, or a matrix with nc columns of counts or proportions. 
#' For either "binomial" or "multinomial", if y is presented as a vector, it will 
#' be coerced into a factor. For family="cox", y should be a two-column matrix with 
#' columns named 'time' and 'status'. The latter is a binary variable, with '1' 
#' indicating death, and '0' indicating right censored. The function Surv() in 
#' package survival produces such a matrix. For family="mgaussian", y is a matrix 
#' of quantitative responses.
#' @param family Response type (see above)
#' @param cores Number of cores used to compute the knockoff statistics by running cv.glmnet.
#' Unless otherwise specified, the number of cores is set equal to two (if available).

#' @param ... additional arguments specific to 'cv.glmnet' (see Details)
#' @return A vector of statistics \eqn{W} (length p)
#' 
#' @details This function uses the \code{glmnet} package to fit a generalized linear model
#' via penalized maximum likelihood.
#' 
#' The knockoff statistics \eqn{W_j} are constructed by taking the difference 
#' between the coefficient of the j-th variable and its knockoff.
#'  
#' By default, the value of the regularization parameter is chosen by 10-fold cross-validation.
#' 
#' The default response family is 'gaussian', for a linear regression model.
#' Different response families (e.g. 'binomial') can be specified by passing an
#' optional parameter 'family'.
#' 
#' The optional \code{nlambda} parameter can be used to control the granularity of the 
#' grid of \eqn{\lambda}'s. The default value of \code{nlambda} is \code{100},
#' where \code{p} is the number of columns of \code{X}.
#' 
#' If the family is 'binomial' and a lambda sequence is not provided by the user, 
#' this function generates it on a log-linear scale before calling 'glmnet'.
#' 
#' For a complete list of the available additional arguments, see \link[glmnet]{cv.glmnet}
#' and \link[glmnet]{glmnet}.
#' 
#' @family statistics for knockoffs
#' 
#' @examples
#' p=100; n=200; k=15
#' mu = rep(0,p); Sigma = diag(p)
#' X = matrix(rnorm(n*p),n)
#' nonzero = sample(p, k)
#' beta = 3.5 * (1:p %in% nonzero)
#' y = X %*% beta + rnorm(n)
#' 
#' knockoffs = function(X) MFKnockoffs.create.gaussian(X, mu, Sigma)
#' # Basic usage with default arguments
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, 
#'                            statistic=MFKnockoffs.stat.glmnet_coef_difference)
#' print(result$selected)
#' 
#' # Advanced usage with custom arguments
#' foo = MFKnockoffs.stat.glmnet_coef_difference
#' k_stat = function(X, X_k, y) foo(X, X_k, y, nlambda=200)
#' result = MFKnockoffs.filter(X, y, knockoffs=knockoffs, statistic=k_stat)
#' print(result$selected)
#' 
#' @rdname MFKnockoffs.stat.lasso_coeff_diff
#' @export
MFKnockoffs.stat.glmnet_coef_difference <- function(X, X_k, y, family='gaussian', cores=2, ...) {
  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  parallel=T
  if (!requireNamespace('doMC', quietly=T)) {
    warning('doMC is not installed. Without parallelization knockoff statistics will be slower to compute', call.=F,immediate.=T)
    parallel=F
  }
  if (!requireNamespace('parallel', quietly=T)) {
    warning('parallel is not installed. Without parallelization knockoff statistics will be slower to compute.', call.=F,immediate.=T)
    parallel=F
  }
    
  # Register cores for parallel computation
  if (parallel) {
    ncores = parallel::detectCores(all.tests = TRUE, logical = TRUE)
    if( cores==2 ) {
      cores = min(2,ncores)
    }
    else {
      if (cores > ncores ) {
        warning(paste("The requested number of cores is not available. Using instead",ncores,"cores"),immediate.=T)
        cores = ncores
      }
    }
    if (cores>1) {
      doMC::registerDoMC(cores=cores)
      parallel = TRUE
    }
    else {
      parallel = FALSE
    }
  }
  
  # Compute statistics
  Z = cv_coeffs_glmnet(cbind(X, X_k), y, family=family, parallel=parallel, ...)
  p = ncol(X)
  orig = 1:p
  abs(Z[orig]) - abs(Z[orig+p])
}

#' @keywords internal
cv_coeffs_glmnet <- function(X, y, nlambda=100, intercept=T, parallel=T, ...) {
  # Standardize variables
  X = normc(X)
  
  n = nrow(X); p = ncol(X)
  n = nrow(X); p = ncol(X)
  
  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family
  
  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }
  
  cv.glmnet.fit <- glmnet::cv.glmnet(X, y, lambda=lambda, intercept=intercept,
                                     standardize=F,standardize.response=F, parallel=parallel, ...)
  
  coef(cv.glmnet.fit, s = "lambda.min")[2:(p+1)]
}
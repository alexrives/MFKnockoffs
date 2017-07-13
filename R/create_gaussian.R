#' Sample multivariate Gaussian knockoff variables
#'
#' Samples multivariate Gaussian fixed-design knockoff variables for the original variables.
#'
#' @param X normalized n-by-p realization of the design matrix
#' @param mu mean vector of length p for X
#' @param Sigma p-by-p covariance matrix for X
#' @param method either 'equi', 'sdp' or 'asdp' (default:'sdp')
#' @param diag_s pre-computed vector of covariances between the original variables and the knockoffs.
#' This will be computed according to 'method', if not supplied
#' @return n-by-p matrix of knockoff variables
#'
#' @family methods for creating knockoffs
#'
#'
#' @references
#'   Candes et al., Panning for Gold: Model-free Knockoffs for High-dimensional Controlled Variable Selection,
#'   arXiv:1610.02351 (2016).
#'   \href{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}{https://statweb.stanford.edu/~candes/MF_Knockoffs/index.html}
#'
#' @export
MFKnockoffs.create.gaussian <- function(X, mu, Sigma, method=c("sdp","asdp","equi"), diag_s=NULL, max_iters = 2500) {
  # Check if covariance matrix if positive-definite
  if (!matrixcalc::is.positive.definite(Sigma)) {
    stop("A positive-definite covariance matrix is required.")
  }
  method = match.arg(method)
  if (is.null(diag_s)) {
    diag_s = diag(switch(match.arg(method),
                    'equi' = MFKnockoffs.knocks.solve_equi(Sigma),
                    'sdp'  = MFKnockoffs.knocks.solve_sdp(Sigma, max_iters = max_iters),
                    'asdp' = MFKnockoffs.knocks.solve_asdp(Sigma, max_iters = max_iters)))
  }
  if (is.null(dim(diag_s))) {
    diag_s = diag(diag_s)
  }

  SigmaInv_s = solve(Sigma,diag_s)
  mu_k = X - sweep(X,2,mu,"-") %*% SigmaInv_s
  Sigma_k = 2*diag_s - diag_s %*% SigmaInv_s
  X_k = mu_k + matrix(rnorm(ncol(X)*nrow(X)),nrow(X)) %*% chol(Sigma_k)
}

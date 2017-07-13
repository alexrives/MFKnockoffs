#' Optimization for SDP knockoffs
#'
#' Solves the optimization problem needed to create SDP knockoffs
#'
#' @param Sigma A positive-definite correlation matrix
#' @param eps Numeric	convergence tolerance for the conic solver (default: 1e-4)
#' @param max_iters The maximum number of iterations for the conic solver (default: 2500)
#' @return The solution \eqn{s} to the semidefinite programming problem defined above
#'
#' @details
#' Solves the semidefinite programming problem:
#'
#'   \deqn{ \mathrm{maximize}      \; \mathrm{sum}(s) \quad
#'           \mathrm{subject} \; \mathrm{to}    0 <= s <= 1, \;
#'                                  2\Sigma - \mathrm{diag}(s) >= 0}
#'
#' If the matrix Sigma supplied by the user is a non-scaled covariance matrix
#' (i.e. its diagonal entries are not all equal to 1), then the appropriate scaling is applied before
#' solving the SDP defined above. The result is then scaled back before being returned, as to match
#' the original scaling of the covariance matrix supplied by the user.
#'
#' @family Optimize knockoffs
#'
#' @export
MFKnockoffs.knocks.solve_sdp <- function(Sigma, eps = 1e-4, max_iters = 2500) {
  if (!requireNamespace('scs', quietly=T))
    stop('scs is not installed', call.=F)
  if (!requireNamespace('Matrix', quietly=T))
    stop('Matrix is not installed', call.=F)

  # Convert the covariance matrix to a correlation matrix
  G = cov2cor(Sigma)
  p = dim(G)[1]

  # Convert problem for SCS

  # First orthant cone constraint
  b1 = rep(0,p)
  A1 = -Matrix::Diagonal(p)

  # Second orthant cone constraint
  b2 = rep(1,p)
  A2 = Matrix::Diagonal(p)

  # Positive-definite cone constraint
  b3 = MFKnockoffs.knocks.vectorize_matrix(2*G)
  d_A3 = MFKnockoffs.knocks.vectorize_matrix(diag(p))
  A3 = Matrix::Diagonal(length(d_A3), x=d_A3)
  A3 = A3[,Matrix::colSums(A3)>0]
  obj = rep(-1,p)

  # Assemble problem
  A = rbind(A1,A2,A3)
  b = c(b1,b2,b3)
  cone = list(l=p*2,s=p)

  # Call SCS with increasing number of iterations until it converges
  converged = F
  converged_tried = 1
  max_iters_run = max_iters
  while (! converged ) {
    control <- list(eps = eps, max_iters = max_iters_run, verbose=F, normalize=F)
    sol <- scs::scs(A, b, obj, cone, control)
    if ( sol$info$statusVal == 1 ) {
      converged = T
    }
    else {
      max_iters_run = max_iters_run * 2
      warning(paste('Conic solver (scs) for SDP knockoffs did not converge after', max_iters_run,
                    'iterations. Trying again with double number of iterations'),immediate.=T)
    }
  }

  # Clip solution to correct numerical errors (domain)
  s = sol$x
  s[s<=0]=1e-8 # Having a zero here would cause instability
  s[s>1]=1

  # Compensate for numerical errors (feasibility)
  psd = 0;
  s_eps = 1e-8;
  while (psd==0) {
    diag_s = diag(s*(1-s_eps))
    GInv_s = solve(G,diag_s)
    Sigma_k = round(2*diag_s - diag_s %*% GInv_s,10)
    if (matrixcalc::is.positive.definite(Sigma_k)) {
      psd  = 1
    }
    else {
      s_eps = s_eps*10
    }
  }
  s = s*(1-s_eps)

  # Scale back the results for a covariance matrix
  return(s*diag(Sigma))
}

#' Vectorize a matrix into the SCS format
#'
#' @rdname vectorize_matrix
#' @keywords internal
MFKnockoffs.knocks.vectorize_matrix = function(M) {
  # Scale the off-diagonal entries by sqrt(2)
  vectorized_matrix = M
  vectorized_matrix[lower.tri(M,diag=FALSE)] = M[lower.tri(M,diag=FALSE)] * sqrt(2)
  # Stack the lower triangular elements column-wise
  vectorized_matrix = vectorized_matrix[lower.tri(vectorized_matrix,diag=TRUE)]
}

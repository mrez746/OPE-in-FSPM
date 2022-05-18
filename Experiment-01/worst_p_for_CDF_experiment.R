worst_p_for_CDF_experiment <- function(pi_l, pi_t, R) {
  # Function for calculating the worst distribution of A*
  # for the case of the CDF estimator
  # Inputs
  #       pi_l (dbl vec) The logging policy
  #       pi_t (dbl vec) The target policy
  #       R (matrix) The reward matrix
  #       n (int) The number of sampled data
  # Outputs:
  #       (dbl vec) the worst p star
  require(CVXR)
  k <- length(pi_t)
  
  compute_A <- function(a) {
    # Function for calculating A_ii entries
    # Inputs:
    #       a (int) Diagonal entry index
    # Output:
    #       (dbl) Value of A_ii
    
    for (i in 1:k) {
      out <- 0
      
      for (j in 1:k) {
        out <- out + R[a, i] * R[a, j] * pi_t[i] * pi_t[j]
      }
    }
    return(out)
  }
  A           <- map_dbl(1:k, ~compute_A(.x))
  A_lag       <- c(0, A[1:(k-1)])
  e_mat       <- matrix(1, k, k)
  indx        <- upper.tri(e_mat)
  e_mat[indx] <- 0
  p           <- Variable(k)
  
  obj <- Maximize(
    sum(
      ((e_mat %*% p - (e_mat %*% p)^2) / pi_l) * A +
        ((e_mat %*% p - (e_mat %*% p)^2) / pi_l) * A_lag
    )
  )
  constr1 <- p >= 0
  constr2 <- sum(p) == 1
  prob <- Problem(obj, constraints = list(constr1, constr2))
  out <- solve(prob)
  out <- out$getValue(p)
  out <- as.vector(out)
  return(out)
}



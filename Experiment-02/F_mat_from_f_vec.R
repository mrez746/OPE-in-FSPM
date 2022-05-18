F_mat_from_f_vec <- function(f1, f0) {
  # make F matrix from f1 and f0 vectors
  # input:
  #       f1 (dbl vec) size k
  #       f0 (dbl vec) size k
  # output:
  #       (matrix) k by k
  k <- length(f1)
  f1_mat <- matrix(rep(f1, k), nrow = k)
  f0_mat <- matrix(rep(f0, k), nrow = k)
  f1_mat[upper.tri(f1_mat, diag = F)] <- 0
  f0_mat[lower.tri(f0_mat, diag = T)] <- 0
  out <- f1_mat + f0_mat
  return(out)
}
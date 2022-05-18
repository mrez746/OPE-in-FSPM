worst_p_star <- function(pi_l, pi_t, F_mat = NULL, imp_sam = F) {
  # function for calculating the worst p_star
  # input:
  #     pi_l (dbl vec) the logging policy (vector of probabilities)
  #     pi_t (dbl vec) the target policy
  #     F_mat (dbl mat) the F matrix
  #     imp_sam (lgl) indicator. TRUE if importance sampling 
  # output:
  #     (dbl vec) the worst p star
  require(CVXR)
  k <- length(pi_l)
  T_mat <- lower.tri(diag(k), diag = T)
  
  if(imp_sam) {
    F_mat <- T_mat * (pi_t / pi_l)
  }
  p <- Variable(k)
  obj <- Minimize((t(pi_t) %*% T_mat %*% p)^2 - t(pi_l) %*% F_mat^2 %*% p)
  constr1 <- p >= 0
  constr2 <- sum(p) == 1
  prob <- Problem(obj, constraints = list(constr1, constr2))
  out <- solve(prob)
  out <- out$getValue(p)
  out <- as.vector(out)
  return(out)
}

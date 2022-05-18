worst_p_star_experiment <- function(pi_l, pi_t, F_mat = NULL, R, imp_sam = F) {
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
  p <- Variable(k)
  obj <- Minimize((t(p) %*% R %*% pi_t)^2 - t(pi_l) %*% F_mat^2 %*% p)
  constr1 <- p >= 0
  constr2 <- sum(p) == 1
  prob <- Problem(obj, constraints = list(constr1, constr2))
  out <- solve(prob)
  out <- out$getValue(p)
  out <- as.vector(out)
  return(out)
}

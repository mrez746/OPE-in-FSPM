calculate_variance <- function(pi_l, pi_t, p_astar, F_mat = NULL, imp_sam = F) {
  # function for calculating the f variance 
  # input:
  #     pi_l (dbl vec) the logging policy (vector of probabilities)
  #     pi_t (dbl vec) the target policy
  #     p_astar (dbl vec) the distribution of the critical value
  #     F_mat (dbl mat) 
  # output:
  #     (dbl) the variance 
  k <- length(pi_l)
  T_mat <- lower.tri(diag(k), diag = T)
  if(imp_sam) {
    F_mat <- T_mat * (pi_t / pi_l)
  }
  out <- t(pi_l) %*% F_mat^2 %*% p_astar - (t(pi_t) %*% T_mat %*% p_astar)^2
  return(out)
}

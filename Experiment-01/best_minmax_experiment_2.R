best_minmax_experiment_2 <- function(pi_l, pi_t, R, ...) {
  # function for finding the best f matrix wrt minmax
  # input:
  #     pi_l (dbl vec) the logging policy (vector of probabilities)
  #     pi_t (dbl vec) the target policy
  #     k (int) cardinality of the action set
  #     ... arguments to pass to CVXR::psolve()
  # output:
  #     (dbl mat) the best f mat
  require(CVXR)
  k <- length(pi_l)
  stopifnot(k == 14)
  # T_mat <- lower.tri(diag(k), diag = T)
  e_vec <- rep(1, k)
  
  A <- matrix(
    rep(pi_l, 2*k), 
    nrow = k, 
    byrow = T
  )
  
  indx_mat <- map(
    0:(k-1), 
    ~tibble(r = .x + 1, c = (1:k) + .x)
  ) %>% 
    bind_rows()
  
  A[as.matrix(indx_mat)] <- 0
  
  alpha <- Variable(1)
  mu <- Variable(1)
  f <- Variable(2*k)
  obj <- Minimize(alpha^2 + mu)
  constr1 <- 2 * alpha %*% R %*% pi_t + mu * e_vec >= A %*% square(f)
  constr2 <- A %*% f == R %*% pi_t
  constr3 <- f[k] == 0
  
  prob <- Problem(
    obj, 
    constraints = list(constr1, constr2, constr3)
  )
  
  out <- psolve(prob, feastol = 1e-3, ...)
  out <- out$getValue(f)
  return(out)
}
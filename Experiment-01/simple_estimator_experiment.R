simple_estimator_experiment <- function(pi_l, pi_t, R, ...) {
  # function for finding the f matrix when f(a,0) is set to 0.
  # Used for the experiment study using real-world data
  # input:
  #     pi_l (dbl vec) the logging policy (vector of probabilities)
  #     pi_t (dbl vec) the target policy
  #     k (int) cardinality of the action set
  #     R (mat) a k_B by k_A matrix containing the reward: R_{a',a} r(a,a')
  #     ... arguments to pass to CVXR::psolve()
  # output:
  #     (list) components:
  #       $f (dbl vec) optimized weights
  #       $t (dbl) from the objective. Here, max regret
  # require(expm)
  # require(CVXR)
  k      <- length(pi_l)
  e_vec  <- rep(1, k-1)
  Pi_inv <- diag(pi_l[-k]^-1)
  
  stopifnot(k == 14)
  
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
  
  out <- solve(A[, (k+1):(2*k)]) %*% R %*% pi_t
  return(out)
  
  # t         <- Variable(0)
  # f         <- Variable(2*k)
  # obj       <- Minimize(0)
  # constr_1 <- A %*% f == R %*% pi_t
  # constr_2 <- f[1] == 0
  # constr_3 <- f[2] == 0
  # constr_4 <- f[3] == 0
  # constr_5 <- f[4] == 0
  # constr_6 <- f[5] == 0
  # constr_7 <- f[6] == 0
  # constr_8 <- f[7] == 0
  # constr_9 <- f[8] == 0
  # constr_10 <- f[9] == 0
  # constr_11 <- f[10] == 0
  # constr_12 <- f[11] == 0
  # constr_13 <- f[12] == 0
  # constr_14 <- f[13] == 0
  # constr_15 <- f[14] == 0
  # 
  # prob <- Problem(
  #   obj, 
  #   constraints = list(
  #     constr_1, constr_2, constr_3, constr_4, constr_5, constr_6,
  #     constr_7, constr_8, constr_9, constr_10, constr_11, constr_12,
  #     constr_13, constr_14, constr_15
  #   )
  # )
  # 
  # rslt <- psolve(prob, feastol = 1e-3, ...)
  # f    <- rslt$getValue(f)
  # t    <- rslt$getValue(t)
  # return(list(f = f, t = t))
}

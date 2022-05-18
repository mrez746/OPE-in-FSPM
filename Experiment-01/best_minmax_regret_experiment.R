best_minmax_regret_experiment <- function(pi_l, pi_t, rho, R, ...) {
  # function for finding the best f matrix wrt minmax
  # for general reward function.
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
  require(expm)
  require(CVXR)
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
  
  B <- rbind(
    diag(k - 1), 
    rep(0, k - 1),
    diag(k - 1),
    -pi_l[-k] / pi_l[k]
  )
  
  D_bar  <- diag(rep(pi_l, 2) * c(k - 1:k, 1:k))
  mat_1  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[1, ])
  mat_2  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[2, ])
  mat_3  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[3, ])
  mat_4  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[4, ])
  mat_5  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[5, ])
  mat_6  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[6, ])
  mat_7  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[7, ])
  mat_8  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[8, ])
  mat_9  <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[9, ])
  mat_10 <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[10, ])
  mat_11 <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[11, ])
  mat_12 <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[12, ])
  mat_13 <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[13, ])
  mat_14 <- sqrtm(Pi_inv - 1) %*% t(B) %*% diag(A[14, ])
  
  t         <- Variable(1)
  f         <- Variable(2*k)
  obj       <- Minimize(t + rho * quad_form(f, D_bar))
  constr_1  <- t >= norm2(mat_1 %*% f)^2
  constr_2  <- t >= norm2(mat_2 %*% f)^2
  constr_3  <- t >= norm2(mat_3 %*% f)^2
  constr_4  <- t >= norm2(mat_4 %*% f)^2
  constr_5  <- t >= norm2(mat_5 %*% f)^2
  constr_6  <- t >= norm2(mat_6 %*% f)^2
  constr_7  <- t >= norm2(mat_7 %*% f)^2
  constr_8  <- t >= norm2(mat_8 %*% f)^2
  constr_9  <- t >= norm2(mat_9 %*% f)^2
  constr_10 <- t >= norm2(mat_10 %*% f)^2
  constr_11 <- t >= norm2(mat_11 %*% f)^2
  constr_12 <- t >= norm2(mat_12 %*% f)^2
  constr_13 <- t >= norm2(mat_13 %*% f)^2
  constr_14 <- t >= norm2(mat_14 %*% f)^2
  constr_15 <- A %*% f == R %*% pi_t
  constr_16 <- f[14] == 0
  
  prob <- Problem(
    obj, 
    constraints = list(
      constr_1, constr_2, constr_3, constr_4, constr_5, constr_6,
      constr_7, constr_8, constr_9, constr_10, constr_11, constr_12,
      constr_13, constr_14, constr_15, constr_16
    )
  )
  
  rslt <- psolve(prob, feastol = 1e-3, ...)
  f    <- rslt$getValue(f)
  t    <- rslt$getValue(t)
  return(list(f = f, t = t))
}

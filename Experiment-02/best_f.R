best_f <- function(pi_l, pi_t, k = 10, ...) {
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
  T_mat <- lower.tri(diag(k), diag = T)
  e_vec <- rep(1, k)
  alpha <- Variable(1)
  mu <- Variable(1)
  F_mat <- Variable(rows = k, cols = k)
  obj <- Minimize(alpha^2 + mu)
  constr1 <- 2 * alpha %*% t(T_mat) %*% pi_t +
    mu * e_vec >= t(square(F_mat)) %*% pi_l
  constr2 <- t(pi_l) %*% F_mat == t(pi_t) %*% T_mat
  
  if (k==10) {
    constr3 <- F_mat[1, 2] == F_mat[1, 3]
    constr4 <- F_mat[1, 3] == F_mat[1, 4]
    constr5 <- F_mat[1, 4] == F_mat[1, 5]
    constr6 <- F_mat[1, 5] == F_mat[1, 6]
    constr7 <- F_mat[1, 6] == F_mat[1, 7]
    constr8 <- F_mat[1, 7] == F_mat[1, 8]
    constr9 <- F_mat[1, 8] == F_mat[1, 9]
    constr10 <- F_mat[1, 9] == F_mat[1, 10]
    
    constr11 <- F_mat[2, 1] == F_mat[2, 2]
    constr12 <- F_mat[2, 3] == F_mat[2, 4]
    constr13 <- F_mat[2, 4] == F_mat[2, 5]
    constr14 <- F_mat[2, 5] == F_mat[2, 6]
    constr15 <- F_mat[2, 6] == F_mat[2, 7]
    constr16 <- F_mat[2, 7] == F_mat[2, 8]
    constr17 <- F_mat[2, 8] == F_mat[2, 9]
    constr18 <- F_mat[2, 9] == F_mat[2, 10]
    
    constr19 <- F_mat[3, 1] == F_mat[3, 2]
    constr20 <- F_mat[3, 2] == F_mat[3, 3]
    constr21 <- F_mat[3, 4] == F_mat[3, 5]
    constr22 <- F_mat[3, 5] == F_mat[3, 6]
    constr23 <- F_mat[3, 6] == F_mat[3, 7]
    constr24 <- F_mat[3, 7] == F_mat[3, 8]
    constr25 <- F_mat[3, 8] == F_mat[3, 9]
    constr26 <- F_mat[3, 9] == F_mat[3, 10]
    
    constr27 <- F_mat[4, 1] == F_mat[4, 2]
    constr28 <- F_mat[4, 2] == F_mat[4, 3]
    constr29 <- F_mat[4, 3] == F_mat[4, 4]
    constr30 <- F_mat[4, 5] == F_mat[4, 6]
    constr31 <- F_mat[4, 6] == F_mat[4, 7]
    constr32 <- F_mat[4, 7] == F_mat[4, 8]
    constr33 <- F_mat[4, 8] == F_mat[4, 9]
    constr34 <- F_mat[4, 9] == F_mat[4, 10]
    
    constr35 <- F_mat[5, 1] == F_mat[5, 2]
    constr36 <- F_mat[5, 2] == F_mat[5, 3]
    constr37 <- F_mat[5, 3] == F_mat[5, 4]
    constr38 <- F_mat[5, 4] == F_mat[5, 5]
    constr39 <- F_mat[5, 6] == F_mat[5, 7]
    constr40 <- F_mat[5, 7] == F_mat[5, 8]
    constr41 <- F_mat[5, 8] == F_mat[5, 9]
    constr42 <- F_mat[5, 9] == F_mat[5, 10]
    
    constr43 <- F_mat[6, 1] == F_mat[6, 2]
    constr44 <- F_mat[6, 2] == F_mat[6, 3]
    constr45 <- F_mat[6, 3] == F_mat[6, 4]
    constr46 <- F_mat[6, 4] == F_mat[6, 5]
    constr47 <- F_mat[6, 5] == F_mat[6, 6]
    constr48 <- F_mat[6, 7] == F_mat[6, 8]
    constr49 <- F_mat[6, 8] == F_mat[6, 9]
    constr50 <- F_mat[6, 9] == F_mat[6, 10]
    
    constr51 <- F_mat[7, 1] == F_mat[7, 2]
    constr52 <- F_mat[7, 2] == F_mat[7, 3]
    constr53 <- F_mat[7, 3] == F_mat[7, 4]
    constr54 <- F_mat[7, 4] == F_mat[7, 5]
    constr55 <- F_mat[7, 5] == F_mat[7, 6]
    constr56 <- F_mat[7, 6] == F_mat[7, 7]
    constr57 <- F_mat[7, 8] == F_mat[7, 9]
    constr58 <- F_mat[7, 9] == F_mat[7, 10]
    
    constr59 <- F_mat[8, 1] == F_mat[8, 2]
    constr60 <- F_mat[8, 2] == F_mat[8, 3]
    constr61 <- F_mat[8, 3] == F_mat[8, 4]
    constr62 <- F_mat[8, 4] == F_mat[8, 5]
    constr63 <- F_mat[8, 5] == F_mat[8, 6]
    constr64 <- F_mat[8, 6] == F_mat[8, 7]
    constr65 <- F_mat[8, 7] == F_mat[8, 8]
    constr66 <- F_mat[8, 9] == F_mat[8, 10]
    
    constr67 <- F_mat[9, 1] == F_mat[9, 2]
    constr68 <- F_mat[9, 2] == F_mat[9, 3]
    constr69 <- F_mat[9, 3] == F_mat[9, 4]
    constr70 <- F_mat[9, 4] == F_mat[9, 5]
    constr71 <- F_mat[9, 5] == F_mat[9, 6]
    constr72 <- F_mat[9, 6] == F_mat[9, 7]
    constr73 <- F_mat[9, 7] == F_mat[9, 8]
    constr74 <- F_mat[9, 8] == F_mat[9, 9]
    
    constr75 <- F_mat[10, 1] == F_mat[10, 2]
    constr76 <- F_mat[10, 2] == F_mat[10, 3]
    constr77 <- F_mat[10, 3] == F_mat[10, 4]
    constr78 <- F_mat[10, 4] == F_mat[10, 5]
    constr79 <- F_mat[10, 5] == F_mat[10, 6]
    constr80 <- F_mat[10, 6] == F_mat[10, 7]
    constr81 <- F_mat[10, 7] == F_mat[10, 8]
    constr82 <- F_mat[10, 8] == F_mat[10, 9]
    constr83 <- F_mat[10, 9] == F_mat[10, 10]
    
    prob <- Problem(
      obj, 
      constraints = list(
        constr1, constr2, constr3, constr4, constr5, constr6,constr7, constr8, 
        constr9, constr10, constr11, constr12, constr13, constr14, constr15, 
        constr16, constr17, constr18, constr19, constr20, constr21, constr22,
        constr23, constr24, constr25, constr26, constr27, constr28, constr29,
        constr30, constr31, constr32, constr33, constr34, constr35, constr36,
        constr37, constr38, constr39, constr40, constr41, constr42, constr43,
        constr44, constr45, constr46, constr47, constr48, constr49, constr50,
        constr51, constr52, constr53, constr54, constr55, constr56, constr57,
        constr58, constr59, constr60, constr61, constr62, constr63, constr64,
        constr65, constr66, constr67, constr68, constr69, constr70, constr71,
        constr72, constr73, constr74, constr75, constr76, constr77, constr78,
        constr79, constr80, constr81, constr82, constr83
      )
    )
  } else {
    stop(paste0("No implementation for k = ", k))
  }
  out <- psolve(prob, feastol = 1e-3, ...)
  out <- out$getValue(F_mat)
  return(out)
}
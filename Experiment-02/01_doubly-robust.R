library(tidyverse)
library(furrr)
source("./worst_p_star.R")
source("./F_mat_from_f_vec.R")
source("./best_minmax_2.R")
source("./best_regret_2.R")
source("./best_minmax_regret.R")

plan(multisession, workers = 4)

df_sim <- readRDS("./clean-data/k10_h1_df_sim.rds")

df_astar_dist <- df_sim %>% 
  filter(id_t == 0) %>% 
  select(id_l, pi_l) %>% 
  rename(
    id_p = id_l, 
    p = pi_l
  )

R_mat <- lower.tri(diag(10), diag = T) * 1
n_vec <- c(50, 1000, 10000)

df_cases_worst <- crossing(
  df_sim,
  tibble(n = n_vec)
)

df_cases_11distr <- crossing(
  df_sim,
  df_astar_dist,
  tibble(n = n_vec)
)

give_rslts <- function(run, add_id_p = FALSE) {
  set.seed(run)
  
  a_vec <- sample(1:10, size = n, prob = pi_l, replace = T)
  b_vec <- sample(1:10, size = n, prob = p, replace = T)
  y_vec <- as.numeric(a_vec >= b_vec)
  
  r_hat <- tapply(y_vec, a_vec, mean)
  r_hat[as.character(setdiff(1:10, names(r_hat)))] <- 0
  stopifnot(length(r_hat) == 10)
  r_hat <- r_hat[as.character(1:10)] %>% as.vector()
  R_dr <- R_mat - r_hat
  
  # f_drmm <- best_f_3(pi_l = pi_l, pi_t = pi_t, R = R_dr)
  # mu_drmm <- f_drmm[y_vec * 10 + a_vec]
  # mu_drmm <- mean(mu_drmm) + as.vector(pi_t %*% r_hat)
  
  mu_mm    <- f_mm[y_vec * 10 + a_vec] %>% mean()
  mu_is    <- f_is[y_vec * 10 + a_vec] %>% mean()
  mu_mmreg <- f_mmreg[y_vec * 10 + a_vec] %>% mean()
  mu_isreg <- f_isreg[y_vec * 10 + a_vec] %>% mean()
  mu_dris  <- f_is[10 + a_vec] * (y_vec - r_hat[a_vec])
  mu_dris  <- mean(mu_dris) + as.vector(pi_t %*% r_hat)
  
  out <- tibble(
    run = run,
    id_l = id_l,
    id_t = id_t,
    n = n,
    mu_is = mu_is,
    mu_mm = mu_mm,
    mu_isreg = mu_isreg,
    mu_mmreg = mu_mmreg,
    mu_dris = mu_dris,
    pi_t_value = pi_t_value
  )
  
  if (add_id_p) out <- out %>% mutate(id_p = id_p)
  
  return(out)
}


for (i in 1:nrow(df_cases_11distr)) {
  p    <- df_cases_11distr$p[[i]]
  pi_l <- df_cases_11distr$pi_l[[i]]
  pi_t <- df_cases_11distr$pi_t[[i]]
  id_l <- df_cases_11distr$id_l[[i]]
  id_t <- df_cases_11distr$id_t[[i]]
  id_p <- df_cases_11distr$id_p[[i]]
  n    <- df_cases_11distr$n[[i]]
  
  f_mm <- best_f_2(pi_l = pi_l, pi_t = pi_t)
  f_is <- c(rep(0, 10), pi_t / pi_l)
  f_mmreg <- best_minmax_regret(pi_l = pi_l, pi_t = pi_t, rho = 0.001)
  f_mmreg <- f_mmreg$f
  f_isreg <- best_regret_2(pi_l = pi_l, pi_t = pi_t)
  
  pi_t_value <- pi_t %*% R_mat %*% p
  pi_t_value <- as.vector(pi_t_value)
  
  out <- future_map(
    1:30,
    ~give_rslts(run = .x, add_id_p = TRUE)
  )
  out <- bind_rows(out)
  print(paste("11 distr i =", i, "of", nrow(df_cases_11distr)))
  
  saveRDS(
    out,
    paste0("Notebooks/12_doubly-robust-estimators/", "dr_11distr_", i, ".rds")
  )
}


for (i in 1:nrow(df_cases_worst)) {
  pi_l <- df_cases_worst$pi_l[[i]]
  pi_t <- df_cases_worst$pi_t[[i]]
  id_l <- df_cases_worst$id_l[[i]]
  id_t <- df_cases_worst$id_t[[i]]
  n    <- df_cases_worst$n[[i]]
  
  f_mm <- best_f_2(pi_l = pi_l, pi_t = pi_t)
  f_is <- c(rep(0, 10), pi_t / pi_l)
  f_mmreg <- best_minmax_regret(pi_l = pi_l, pi_t = pi_t, rho = 0.001)
  f_mmreg <- f_mmreg$f
  f_isreg <- best_regret_2(pi_l = pi_l, pi_t = pi_t)
  F_mm <- F_mat_from_f_vec(f1 = f_mm[11:20], f0 = f_mm[1:10])
  F_mmreg <- F_mat_from_f_vec(f1 = f_mmreg[11:20], f0 = f_mmreg[1:10])
  F_isreg <- F_mat_from_f_vec(f1 = f_isreg[11:20], f0 = f_isreg[1:10])
  
  worst_p_mm <- worst_p_star(pi_l = pi_l, pi_t = pi_t, F_mat = F_mm)
  worst_p_is <- worst_p_star(pi_l = pi_l, pi_t = pi_t, imp_sam = T)
  worst_p_mmreg <- worst_p_star(pi_l = pi_l, pi_t = pi_t, F_mat = F_mmreg)
  worst_p_isreg <- worst_p_star(pi_l = pi_l, pi_t = pi_t, F_mat = F_isreg)
  
  nms_worst <- c("mm", "is", "mmreg", "isreg")
  
  out_2 <- vector("list", length(nms_worst))
  
  for (j in 1:length(nms_worst)) {
    p <- get(paste0("worst_p_", nms_worst[j]))
    p[p < 1e-10] <- 0
    pi_t_value <- pi_t %*% R_mat %*% p
    pi_t_value <- as.vector(pi_t_value)
    
    out <- future_map(
      1:30,
      ~give_rslts(run = .x)
    )
    
    out_2[[j]] <- out %>% 
      bind_rows() %>% 
      mutate(worst_p = nms_worst[j])
  }
  out_2 <- bind_rows(out_2)
  print(paste("Worst i =", i, "of", nrow(df_cases_worst)))
  
  saveRDS(
    out_2,
    paste0("Notebooks/12_doubly-robust-estimators/", "dr_worst_", i, ".rds")
  )
}

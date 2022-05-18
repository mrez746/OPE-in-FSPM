library(tidyverse)
library(furrr)
source("./Notebooks/10_real-world-example/experiment-functions.R")
source("./CDF_est_method.R")
source("./CDF_est_method_2.R")

experiment_list <- readRDS(
  "./Notebooks/10_real-world-example/experiment-policies-weights.rds"
)
weights_df   <- experiment_list$weights_df
tree_target  <- experiment_list$tree_target
tree_logging <- experiment_list$tree_logging
df_test      <- experiment_list$df_test
R            <- experiment_list$R

g_target_test  <- find_leaf_of_point(tree_target, df_test)
g_logging_test <- find_leaf_of_point(tree_logging, df_test)

p_astar_df <- readRDS(
  "./Notebooks/10_real-world-example/experiment_worst_p_all.rds"
)

p_astar_df <- p_astar_df %>% 
  mutate(
    p_astar = map(p_astar, ~if_else(.x < 1e-20, 0, .x)), 
    p_astar = map(p_astar, ~.x/sum(.x))
  )

# not all combination of pi_t and pi_l nodes exists
p_astar_df <- tibble(
  g_logging = g_logging_test,
  g_target  = g_target_test
) %>% 
  count(g_logging, g_target) %>% 
  select(-n) %>% 
  left_join(p_astar_df, by = c("g_logging", "g_target"))

node_w <- tibble(
  g_logging = g_logging_test, 
  g_target = g_target_test
) %>% 
  count(g_logging, g_target) %>% 
  mutate(
    node_w = n / sum(n)
  ) %>% 
  select(-n)

run_experiment <- function(seed, worst_p_estimator, 
                           w1, w2, subsample_n = NULL) {
  # Function for comparing the mm, mmregret, simple and CDF method on 
  # a bootstrap data set
  # Inputs:
  #       seed (int) seed for making the bootstrap sample
  #       w1 (dbl) w1*empirical_PMF + (1-w1)*rev_empirical_PMF for pi_l
  #       w2 (dbl) w2*empirical_PMF + (1-w2)*rev_empirical_PMF for pi_t
  #       worst_p_estimator (chr) worst A* from which estimator
  #       subsample_n (int) # data points to subsample 
  # Env. objects:
  #       df_test (dataframe) 
  #       weights_df (dataframe) of policies and weights of estimators. Columns:
  #         $ w_log (dbl)
  #         $ w_tar (dbl) 
  #         $ estimator (chr)  One of: mm, mmregret, simple
  #         $ g_logging (int)  The leaf ID of the tree of pi_l
  #         $ g_target (int)   The leaf ID of the tree of pi_t
  #         $ f0 (list)        14-dimensional vector of calculated weights
  #         $ f1 (list)        14-dimensional vector of calculated weights
  #         $ p_logging (list) 14-dimensional vectors of probs.
  #         $ p_target (list)  14-dimensional vectors of probs.
  #       p_astar_df (dataframe) of worst A* PMFs
  #       node_w (dataframe) node weights from the tree of pi_t
  #       find_leaf_of_point (fun) 
  # Outputs:
  #       (dataframe) w/ components:
  #         $ seed (int)
  #         $ pi_t_val (dbl)     The true pi_t value
  #         $ est_mm (dbl)       Estimate of pi_t value by mm
  #         $ est_mmregret (dbl) Estimate of pi_t value by mmregret
  #         $ est_simple (dbl)   Estimate of pi_t value by simple
  #         $ est_CDF (dbl)      Estimate of pi_t value by CDF
  set.seed(seed)
  m <- nrow(df_test)
  
  if (is.null(subsample_n)) {
    boot_i <- sample(1:m, m, replace = T) # bootstrap
  } else {
    boot_i <- sample(1:m, subsample_n, replace = F)
  }
  
  df <- tibble(
    g_logging = g_logging_test[boot_i],
    g_target  = g_target_test[boot_i]
  )
  
  this_p_astar_df <- p_astar_df %>% 
    filter(w_log == w1, w_tar == w2, estimator == worst_p_estimator) %>% 
    select(g_logging, g_target, p_astar)
  
  weights_mm <- weights_df %>%
    filter(w_log == w1, w_tar == w2, estimator == "mm") %>%
    select(g_logging, g_target, f0, f1)
  
  weights_mmregret <- weights_df %>%
    filter(w_log == w1, w_tar == w2, estimator == "mmregret") %>%
    select(g_logging, g_target, f0, f1)
  
  weights_simple <- weights_df %>%
    filter(w_log == w1, w_tar == w2, estimator == "simple") %>%
    select(g_logging, g_target, f0, f1)
  
  pi_l <- weights_df %>% 
    filter(w_log == w1, w_tar == w2, estimator == "mm", g_target == 1) %>% 
    select(g_logging, p_logging)
  
  pi_t <- weights_df %>% 
    filter(w_log == w1, w_tar == w2, estimator == "mm", g_logging == 1) %>% 
    select(g_target, p_target)
  
  pi_t_val <- pi_t %>% 
    left_join(this_p_astar_df, by = "g_target") %>% 
    mutate(
      node_val = map2_dbl(p_astar, p_target, ~t(.x) %*% R %*% .y)
    ) %>% 
    left_join(node_w, by = c("g_target", "g_logging")) %>% 
    mutate(
      node_est_weighted = node_w * node_val
    )
  pi_t_val <- sum(pi_t_val$node_est_weighted)
  
  df <- df %>% 
    left_join(pi_l, by = "g_logging") %>% 
    left_join(this_p_astar_df, by = c("g_logging", "g_target")) %>% 
    mutate(
      a = map_int(p_logging, ~sample(1:14, 1, prob = .x)),
      astar = map_int(p_astar, ~sample(1:14, 1, prob = .x)),
      obs = as.numeric(a >= astar)
    )
  
  est_mm <- df %>%
    left_join(weights_mm, by = c("g_logging", "g_target")) %>%
    mutate(
      est_0 = map2_dbl(a, f0, ~.y[.x]),
      est_1 = map2_dbl(a, f1, ~.y[.x]),
      single_est = if_else(obs == 0, est_0, est_1)
    )
  est_mm <- mean(est_mm$single_est)
  
  est_mmregret <- df %>%
    left_join(weights_mmregret, by = c("g_logging", "g_target")) %>%
    mutate(
      est_0 = map2_dbl(a, f0, ~.y[.x]),
      est_1 = map2_dbl(a, f1, ~.y[.x]),
      single_est = if_else(obs == 0, est_0, est_1)
    )
  est_mmregret <- mean(est_mmregret$single_est)
  
  est_simple <- df %>%
    left_join(weights_simple, by = c("g_logging", "g_target")) %>%
    mutate(
      est_0 = map2_dbl(a, f0, ~.y[.x]),
      est_1 = map2_dbl(a, f1, ~.y[.x]),
      single_est = if_else(obs == 0, est_0, est_1)
    )
  est_simple <- mean(est_simple$single_est)
  
  est_CDF <- CDF_est_method(
    sample_df = df,
    pi_t = pi_t %>% mutate(a = list(1:14)) %>% unnest
  )
  
  est_CDF_2 <- CDF_est_method_2(
    sample_df = df, 
    pi_t = pi_t,
    pi_l = pi_l
  )
  
  out <- tibble(seed = seed,
                pi_t_val = pi_t_val,
                est_mm = est_mm,
                est_mmregret = est_mmregret,
                est_simple = est_simple,
                est_CDF_reg = est_CDF,
                est_CDF_IS = est_CDF_2)
  return(out)
}


df_par <- crossing(
  subsample_n = c(1000, 100, 50),
  estimator = c("mm", "mmregret", "CDF"),
  w2 = seq(0, 1, by = .2)
)

plan(multisession, workers = 8)

for (i in 1:nrow(df_par)) {
  
  tictoc::tic()
  
  boot_rslt <- future_map(
    1:200,
    ~run_experiment(
      seed = .x, 
      worst_p_estimator = df_par$estimator[i], 
      w1 = 1, 
      w2 = df_par$w2[i], 
      subsample_n = df_par$subsample_n[i]
    )
  )
  
  boot_rslt <- bind_rows(boot_rslt) %>% 
    mutate(
      p_astar_estimator = df_par$estimator[i],
      w2 = df_par$w2[i],
      subsample_n = df_par$subsample_n[i]
    )
  
  sv_nm <- paste0(
    "./Notebooks/10_real-world-example/boot-rslt_batch-3/",
    "boot_rslt_", i, ".rds"
  )
  saveRDS(boot_rslt, sv_nm)
  
  print(paste(
    "--- i", i, "of", nrow(df_par)
  ))
  tictoc::toc()
}




df_par <- crossing(
  subsample_n = c(7000),
  estimator = c("mm", "mmregret", "CDF"),
  w2 = seq(0, 1, by = .2)
)
plan(multisession, workers = 8)
seed_vec <- seq(1, 200 * 10 + 1, by = 200)

for (j in 1:length(seed_vec)) {
  seeds <- seed_vec[j]:(seed_vec[j] + 199)
  out <- vector("list", nrow(df_par))
  tictoc::tic()
  
  for (i in 1:nrow(df_par)) {
    
    boot_rslt <- future_map(
      seeds,
      ~run_experiment(
        seed = .x, 
        worst_p_estimator = df_par$estimator[i], 
        w1 = 1, 
        w2 = df_par$w2[i], 
        subsample_n = df_par$subsample_n[i]
      )
    )
    
    boot_rslt <- bind_rows(boot_rslt) %>% 
      mutate(
        p_astar_estimator = df_par$estimator[i],
        w2 = df_par$w2[i],
        subsample_n = df_par$subsample_n[i],
        replication = j
      )
    out[[i]] <- boot_rslt
  }
  out <- bind_rows(out)
  
  print(paste(
    "--- j", j, "of", length(seed_vec)
  ))
  tictoc::toc()
  
  sv_nm <- paste0(
    "./Notebooks/10_real-world-example/boot-rslt_batch-4/",
    "boot_rslt_", j, ".rds"
  )
  saveRDS(out, sv_nm)
}






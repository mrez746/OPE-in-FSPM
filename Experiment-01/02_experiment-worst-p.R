library(tree)
library(Matrix)
library(CVXR)
library(tidyverse)

source("./F_mat_from_f_vec.R")
source("./worst_p_star_experiment.R")
source("./worst_p_for_CDF_experiment.R")

experiment_list <- readRDS(
  "./Notebooks/10_real-world-example/experiment-policies-weights.rds"
)
df <- experiment_list$weights_df
R  <- experiment_list$R

df <- df %>% 
  filter(estimator %in% c("mm", "mmregret")) %>% 
  mutate(
    F_mat = map2(f1, f0, ~F_mat_from_f_vec(f1 = .x, f0 = .y))
  )
out <- vector("list", nrow(df))

for (i in 1:length(out)) {
  pi_l  <- df$p_logging[[i]]
  pi_t  <- df$p_target[[i]]
  F_mat <- df$F_mat[[i]]
  
  out[[i]] <- worst_p_star_experiment(
    pi_l = pi_l, 
    pi_t = pi_t, 
    F_mat = F_mat, 
    R = R
  )
  print(paste("i", i, "of", length(out)))
}

out <- map2(
  out, 
  df$id, 
  ~tibble(worst_p = list(.x),
          id = .y)
)
out <- bind_rows(out)

# saveRDS(
#   out,
#   "./Notebooks/10_real-world-example/experiment-worst_p.rds"
# )

df2 <- df %>% 
  group_by(g_logging, g_target, w_log, w_tar) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(id, g_logging, p_logging, g_target, p_target)
# n <- experiment_list$df_test %>% nrow()

df2 <- df2 %>% 
  mutate(
    worst_p_CDF = map2(
      p_logging,
      p_target,
      ~worst_p_for_CDF_experiment(
        pi_l = .x, 
        pi_t = .y, 
        R = R
      )
    )
  )

# saveRDS(
#   df2,
#   "./Notebooks/10_real-world-example/experiment_worst_p_CDF.rds"
# )

df1 <- readRDS(
  "./Notebooks/10_real-world-example/experiment-worst_p.rds"
)

df_cdf <- readRDS(
  "./Notebooks/10_real-world-example/experiment_worst_p_CDF.rds"
)


df_cdf <- df_cdf %>% 
  left_join(
    df %>% select(id, w_log, w_tar), 
    by = c("id")
  ) %>%
  select(-id, -p_logging, -p_target) %>% 
  rename(p_astar = worst_p_CDF) %>% 
  mutate(estimator = "CDF")

df1 <- df1 %>% 
  left_join(df, by = "id")

df1 <- df1 %>% 
  rename(p_astar = worst_p) %>% 
  select(g_logging, g_target, w_log, w_tar, p_astar, estimator)

worst_p_df <- bind_rows(df1, df_cdf)

# saveRDS(
#   worst_p_df,
#   "./Notebooks/10_real-world-example/experiment_worst_p_all.rds"
# )


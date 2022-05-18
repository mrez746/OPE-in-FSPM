library(tree)
library(Matrix)
library(CVXR)
library(tidyverse)

source("./best_minmax_regret_experiment.R")
source("./best_minmax_experiment_2.R")
source("./simple_estimator_experiment.R")
source("./CDF_est_method.R")
source("./Notebooks/10_real-world-example/experiment-functions.R")

test_percent <- 70

df <- read_csv("./data/dataset_diabetes/diabetic_data.csv")
names(df) <- str_replace_all(names(df), "-", "_")

circulatory     <- c(390:459, 785) %>% as.character()
respiratory     <- c(460:519, 786) %>% as.character()
digestive       <- c(520:579, 787) %>% as.character()
injury          <- 800:999 %>% as.character()
musculoskeletal <- 710:739 %>% as.character()
genitourinary   <- c(580:629, 788) %>% as.character()
neoplasms       <- 140:239 %>% as.character()

diabetes <- df$diag_1 %>% 
  unique %>% 
  str_subset("250\\.")

df <- df %>% 
  mutate(
    diag_group = case_when(
      diag_1 %in% diabetes ~ "diabetes",
      diag_1 %in% circulatory ~ "circulatory",
      diag_1 %in% respiratory ~ "respiratory",
      diag_1 %in% digestive ~ "digestive",
      diag_1 %in% injury ~ "injury",
      diag_1 %in% musculoskeletal ~ "musculoskeletal",
      diag_1 %in% genitourinary ~ "genitourinary",
      diag_1 %in% neoplasms ~ "neoplasms",
      TRUE ~ "other"
    )
  )

df_sml <- df %>% 
  select(
    time_in_hospital,
    diag_group, A1Cresult, age#, medical_specialty
  ) %>% 
  mutate(
    diag_group = as.factor(diag_group),
    A1Cresult = as.factor(A1Cresult),
    age = as.factor(age),
    #medical_specialty = as.factor(medical_specialty)
  )
R <- matrix(0, nrow = 14, ncol = 14)

for (i in 1:14) {
  for (j in 1:14) {
    R[i, j] = -abs(i - j)
  }
}

set.seed(999)
n        <- nrow(df_sml)
train_i  <- sample(1:n, floor(n * (100 - test_percent) / 100))
df_train <- df_sml[train_i, ]
df_test  <- df_sml[-train_i, ]

tree_logging <- tree::tree(
  time_in_hospital ~ ., 
  data = df_train, #%>% select(-medical_specialty),
  control = tree::tree.control(
    nrow(df_train), 
    mincut = 5, 
    minsize = 10,
    mindev = 0.003)
)

tree_target <- tree::tree(
  time_in_hospital ~ ., 
  data = df_train, #%>% select(-medical_specialty),
  control = tree::tree.control(
    nrow(df_train), 
    mincut = 5, 
    minsize = 10,
    mindev = 0.001)
)

policy_logging_good <- make_policy(
  tree_logging, 
  df_train
)

policy_target_good <- make_policy(
  tree_target, 
  df_train
)

policy_logging_bad <- policy_logging_good %>% 
  group_by(g) %>% 
  mutate(p = 1/p, p = p / sum(p)) %>% 
  ungroup()

policy_target_bad <- policy_target_good %>% 
  group_by(g) %>% 
  mutate(p = 1/p, p = p / sum(p)) %>% 
  ungroup()

w1 <- seq(0, 1, by = .2)
w2 <- seq(0, 1, by = .2)
out1 <- vector("list", length(w1))

for (i in 1:length(w1)) {
  
  pi_l <- left_join(
    policy_logging_good %>% select(-n) %>% rename(p_good = p),
    policy_logging_bad %>% select(-n) %>% rename(p_bad = p),
    by = c("g", "time_in_hospital")
  ) %>% 
    mutate(p = w1[i] * p_good + (1 - w1[i]) * p_bad)
  
  logging_p_df <- pi_l %>% 
    group_by(g) %>% 
    arrange(time_in_hospital) %>% 
    nest() %>% 
    mutate(p_logging = map(data, "p")) %>% 
    rename(g_logging = g) %>% 
    select(-data)
  
  out2 <- vector("list", length(w2))
  
  for (j in 1:length(w2)) {
    
    pi_t <- left_join(
      policy_target_good %>% select(-n) %>% rename(p_good = p),
      policy_target_bad %>% select(-n) %>% rename(p_bad = p),
      by = c("g", "time_in_hospital")
    ) %>% 
      mutate(
        p = w2[j] * p_good + (1 - w2[j]) * p_bad
      )
    
    target_p_df <- pi_t %>% 
      group_by(g) %>% 
      arrange(time_in_hospital) %>% 
      nest() %>% 
      mutate(p_target = map(data, "p")) %>% 
      rename(g_target = g) %>% 
      select(-data)
    
    weights_df <- crossing(
      logging_p_df,
      target_p_df
    )
    
    weights_df_mmregret <- weights_df %>%
      mutate(
        optimization_rslt = map2(
          p_logging,
          p_target,
          ~best_minmax_regret_experiment(
            pi_l = .x,
            pi_t = .y,
            rho = 0.001,
            R = R
          )
        )
      )
    
    weights_df_mmregret <- weights_df_mmregret %>%
      mutate(
        f0 = map(optimization_rslt, ~.x$f[1:14, 1]),
        f1 = map(optimization_rslt, ~.x$f[15:28, 1]),
        estimator = "mmregret"
      )
    
    weights_df_mm <- weights_df %>%
      mutate(
        optimization_rslt = map2(
          p_logging,
          p_target,
          ~best_minmax_experiment_2(
            pi_l = .x,
            pi_t = .y,
            R = R
          )
        )
      )
    
    weights_df_mm <- weights_df_mm %>%
      mutate(
        f0 = map(optimization_rslt, ~.x[1:14, 1]),
        f1 = map(optimization_rslt, ~.x[15:28, 1]),
        estimator = "mm"
      )
    
    weights_df_simple <- weights_df %>%
      mutate(
        optimization_rslt = map2(
          p_logging,
          p_target,
          ~simple_estimator_experiment(
            pi_l = .x,
            pi_t = .y,
            R = R
          )
        )
      )
    
    weights_df_simple <- weights_df_simple %>%
      mutate(
        f0 = list(rep(0, 14)),
        f1 = map(optimization_rslt, ~as.vector(.x)),
        estimator = "simple"
      )
    
    out2[[j]] <- weights_df_mmregret %>% 
      bind_rows(weights_df_mm) %>% 
      bind_rows(weights_df_simple) %>% 
      mutate(
        w_log = w1[i],
        w_tar = w2[j]
      )
    print(paste(
      "i", i, "of", length(w1), "---",
      "j", j, "of", length(w2)
    ))
  }
  out1[[i]] <- bind_rows(out2)
}
out1 <- bind_rows(out1) %>% 
  mutate(id = row_number()) %>% 
  select(id, everything())

out <- list(
  weights_df = out1,
  tree_logging = tree_logging, 
  tree_target = tree_target,
  policy_logging_good = policy_logging_good, 
  policy_target_good = policy_target_good,
  policy_logging_bad = policy_logging_bad, 
  policy_target_bad = policy_target_bad,
  df_test = df_test,
  R = R, 
  df_train = df_train
)

saveRDS(
  out,
  "./Notebooks/10_real-world-example/experiment-policies-weights.rds"
)

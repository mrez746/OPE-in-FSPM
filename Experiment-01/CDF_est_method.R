CDF_est_method <- function(sample_df, pi_t) {
  # function for estimating the target policy value by estimating the 
  # CDF of the critical value, for the experimental study
  # assumes k = 14 and input data with certain attributes
  # input
  #     sample_df (df) w/ columns:
  #         $a
  #         $g_target
  #         $obs
  #     pi_t (df) 2/ columns:
  #         $a
  #         $g_target
  #         $p_target
  # output
  #     (dbl) target value estimate
  F_hat <- sample_df %>% 
    group_by(g_target, a) %>% 
    summarise(
      F_a = mean(obs), 
      .groups = "drop"
    )
  
  F_hat_zero <- crossing(
    a = unique(pi_t$a),
    g_target = unique(pi_t$g_target)
  ) %>% 
    anti_join(F_hat, by = c("a", "g_target")) %>% 
    mutate(F_a = 0)
  
  F_hat <- F_hat %>% 
    bind_rows(F_hat_zero) %>% 
    group_by(g_target) %>% 
    arrange(a) %>% 
    mutate(
      p_a = F_a - lag(F_a), 
      p_a = if_else(is.na(p_a), F_a, p_a)
    ) %>% 
    ungroup() %>% 
    arrange(g_target, a)
  
  F_hat <- F_hat %>% 
    rename(
      aprime = a,
      F_aprime = F_a,
      p_aprime = p_a
    )
  
  g_aprime <- pi_t %>% 
    select(g_target, a, p_target)
  
  g_aprime <- crossing(
    g_aprime,
    tibble(aprime = 1:14)
  )
  
  g_aprime <- g_aprime %>% 
    group_by(g_target, aprime) %>% 
    summarise(
      g = sum(-abs(a - aprime) * p_target),
      .groups = "drop"
    )
  
  simple_val_df <- F_hat %>% 
    left_join(
      g_aprime, 
      by = c("g_target", "aprime")
    ) %>% 
    group_by(g_target) %>% 
    summarise(
      val = sum(p_aprime * g),
      .groups = "drop"
    )
  
  # mean(simple_val_df$val)
  
  weight_df <- sample_df %>% 
    count(g_target) %>% 
    mutate(g_weight = n / sum(n))
  
  weight_df_zero <- tibble(
    g_target = unique(pi_t$g_target),
    g_weight = 0
  ) %>% 
    anti_join(weight_df, by = "g_target")
  
  weight_df <- bind_rows(weight_df, weight_df_zero)
  
  simple_val_df <- simple_val_df %>% 
    left_join(
      weight_df, 
      by = "g_target"
    ) %>% 
    mutate(
      weighted_val = val * g_weight
    )
  
  out <- sum(simple_val_df$weighted_val)
  return(out)
}


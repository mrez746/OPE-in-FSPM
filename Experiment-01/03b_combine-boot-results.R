library(tidyverse)

combine_files <- function(dir_nm) {
  
  files <- dir(paste0(
    "./Notebooks/10_real-world-example/", dir_nm
  ))
  out <- vector("list", length(files))
  
  for (i in 1:length(out)) {
    
    out[[i]] <- readRDS(paste0(
      "./Notebooks/10_real-world-example/", 
      dir_nm, "/",
      files[i]
    ))
  }
  out <- bind_rows(out)
  return(out)
}

out1 <- combine_files("boot-rslt_batch-1")
out2 <- combine_files("boot-rslt_batch-2")
out3 <- combine_files("boot-rslt_batch-3")
out4 <- combine_files("boot-rslt_batch-4")
  
out <- out1 %>% 
  mutate(est_CDF_reg = est_CDF) %>% 
  left_join(
    out2 %>% select(seed, p_astar_estimator, w2, est_CDF_2),
    by = c("seed", "p_astar_estimator", "w2")
  ) %>% 
  mutate(est_CDF_IS = est_CDF_2) %>% 
  bind_rows(out3) %>% 
  bind_rows(
    out4 %>% filter(replication == 1) %>% select(-replication)
  ) %>% 
  mutate(
    subsample_n = if_else(is.na(subsample_n),
                          71237,
                          subsample_n)
  )

saveRDS(
  out,
  "./Notebooks/10_real-world-example/boot-rslts.rds"
)  

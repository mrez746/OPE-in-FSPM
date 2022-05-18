library(tidyverse)

nms <- dir("Notebooks/12_doubly-robust-estimators/output-data/") %>% 
  str_subset("dr_.+rds")

out <- vector("list", length(nms))

for (i in 1:length(out)) {
  
  df <- readRDS(paste0("Notebooks/12_doubly-robust-estimators/output-data/", 
                       nms[i]))
  
  if (str_detect(nms[i], "worst")) {
    df <- df %>% rename(id_p = worst_p)
  } else{
    df <- df %>% mutate(id_p = as.character(id_p))
  }
  out[[i]] <- df
  print(paste(i, "of", length(out)))
}
df_comparison <- bind_rows(out)

saveRDS(df_comparison, 
        "Notebooks/12_doubly-robust-estimators/df_comparison.rds")

df_comparison


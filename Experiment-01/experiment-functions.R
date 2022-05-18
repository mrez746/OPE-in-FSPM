find_leaf_of_point <- function(mod_tree, df) {
  # find which leaf points in df fall into
  leaves_vals <- mod_tree$frame %>% 
    filter(var == "<leaf>") %>% 
    .$yval
  
  if (any(duplicated(leaves_vals))) {
    error("Error: Duplicated leave values")
  }
  
  out <- map_int(
    tree:::predict.tree(mod_tree, df), 
    ~which(.x == leaves_vals)
  )
  return(out)
}

make_policy <- function(mod_tree, df) {
  # calculate the conditional pmf using the data points in each leaf
  g <- find_leaf_of_point(mod_tree, df)
  
  out <- df %>% 
    mutate(g = g) %>% 
    count(g, time_in_hospital) %>% 
    group_by(g) %>% 
    mutate(p = n/sum(n)) %>% 
    ungroup()
  return(out)
}
#' Cell engrafment time frames as in https://pmc.ncbi.nlm.nih.gov/articles/PMC6467214/
#' - “initial engraftment” 0–3 months
#' - “stabilization” 3 months to 1 year
#' - “homeostasis” >1 year
#' - “long-term” engraftment >2 years after transplant.
engrafment_phase <- Vectorize(function(tp) {
  
  if(tp < 90) {
    return("Initial Engrafment")
  }
  
  if(tp < 365) {
    return("Stabilization")
  }
  
  if(tp < 2*365) {
    return("Homeostasis")
  }
  
  return("Long-term Engrafment")
})

get_diversity_df <- function(sim_data, tps = NULL) {
  if(is.null(tps)) {
    tps <- as.integer(substring(colnames(sim_data), 2))
  }
  
  nr_cells <- colSums(sim_data[, paste0("t", tps)])
  
  eventually_appear <- sim_data[, paste0("t", tps)] > 0
  for(idx in 1:(length(tps) - 1)) {
    eventually_appear[, length(tps) - idx] <-
      eventually_appear[, length(tps) - idx] |
      eventually_appear[, length(tps) - idx + 1]
  }
  nr_clones <- colSums(eventually_appear)
  
  return(data.frame(
    Time = tps,
    Nr_Cells = nr_cells,
    Nr_Clones = nr_clones,
    Mean_Size = nr_cells/nr_clones
  ))
}

get_contribution_dist_df <- function(sim_data, tps = NULL, break_list = NULL) {
  if(is.null(tps)) {
    tps <- as.integer(substring(colnames(sim_data), 2))
  }
  
  if(is.null(break_list)) {
    break_list <- rep("Sturges", length(tps))
    names(break_list) <- paste0("t", tps)
  }
  
  return(data.frame(do.call(rbind, lapply(tps, function(tp) {
    curr_data <- sim_data[, paste0("t", tp)]
    
    existing_clones <- which(curr_data > 0)
    existing_clone_data <- curr_data[existing_clones]
    
    cell_count <- sum(existing_clone_data)
    hist_data <- hist(
      log10(existing_clone_data/cell_count), 
      plot = FALSE,
      breaks = break_list[[paste0("t", tp)]]
    )
    non_zero <- which(hist_data$counts > 0)
    
    return(data.frame(
      Time = tp,
      Engrafment_Phase = engrafment_phase(tp),
      
      Clone_Contribution = 10**hist_data$breaks[non_zero],
      Frequency = hist_data$counts[non_zero]/length(existing_clones)
    ))
  }))))
}

get_so_df <- function(sim_data, tps = NULL) {
  if(is.null(tps)) {
    tps <- as.integer(sapply(colnames(sim_data), function(t_name) substring(t_name, 2)))
  }
  
  clone_occurrences <- rowSums(sim_data[, paste0("t", tps)] > 0)
  so_clones <- which(clone_occurrences == 1)
  total_cells <- colSums(sim_data[, paste0("t", tps)])
  so_cells <- colSums(sim_data[so_clones, paste0("t", tps)])
  
  nr_clones <- colSums(sim_data[, paste0("t", tps)] > 0)
  nr_so_clones <- colSums(sim_data[so_clones, paste0("t", tps)] > 0)
  
  return(data.frame(
    SO_Cells = so_cells,
    Total_Cells = total_cells,
    SO_Cell_Contribution = so_cells/total_cells,
    SO_Clones = nr_so_clones,
    Total_Clones = nr_clones,
    SO_Clone_Contribution = nr_so_clones/nr_clones,
    Time = tps
  ))
}

get_full_diffs_df <- function(sim_diffs) {
  df <- bind_rows(
    #### Z13264 ####
    Z13264 = bind_rows(
      Diversity = sim_diffs$Z13264$Diversity %>% rename(
        X = Time,
        Y = Nr_Clones_Sim
      ),
      SO = sim_diffs$Z13264$SO %>% rename(
        X = Time,
        Y = SO_Clone_Contribution_Sim
      ),
      Contribution_397 = sim_diffs$Z13264$Contribution %>%
        filter(Time == 397) %>%
        rename(
          X = Clone_Contribution,
          Y = Frequency_Sim
        ),
      Contribution_462 = sim_diffs$Z13264$Contribution %>%
        filter(Time == 462) %>%
        rename(
          X = Clone_Contribution,
          Y = Frequency_Sim
        ),
      .id = "Graph"
    ),
    
    #### Z14004 ####
    Z14004 = bind_rows(
      Diversity = sim_diffs$Z14004$Diversity %>% rename(
        X = Time,
        Y = Nr_Clones_Sim
      ),
      SO = sim_diffs$Z14004$SO %>% rename(
        X = Time,
        Y = SO_Clone_Contribution_Sim
      ),
      Contribution_362 = sim_diffs$Z14004$Contribution %>%
        filter(Time == 362) %>%
        rename(
          X = Clone_Contribution,
          Y = Frequency_Sim
        ),
      Contribution_495 = sim_diffs$Z14004$Contribution %>%
        filter(Time == 495) %>%
        rename(
          X = Clone_Contribution,
          Y = Frequency_Sim
        ),
      
      # End ####
      .id = "Graph"
    ),
    
    .id = "Subject"
  ) %>% select(
    Subject, Graph,
    X, Y, Res
  )
  
  return(df)
}
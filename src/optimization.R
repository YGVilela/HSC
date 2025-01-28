load("data/derived_radtke.Rda")
source("src/utils.R")
source("src/sampling.R")
source("src/Simplified_HSC_Model.R")

library(tidyverse)

# Base optimization params ####
base_params = c(
  pA = 0.8,
  dA = 0.08,
  tQA = 0.02,
  tAQ = 0.02,
  pQ = 0,
  q_proportion = 0.5,
  clone_mult = 3,
  space_mult_A = 1.5,
  space_mult_Q = 2.5
)

# Simulation constants ####

# Get time and sample sizes from data
distribution_time <- list(
  "Z13264" = 249,
  "Z14004" = 362
)

tps <- list(
  "Z13264" = (clones_over_time %>% filter(Subject == "Z13264"))$Time,
  "Z14004" = (clones_over_time %>% filter(Subject == "Z14004"))$Time
)
all_tps <- sort(clones_over_time$Time)
data_sample_sizes <- list(
  "Z13264" = clone_data %>%
    filter(Subject == "Z13264") %>%
    group_by(Time) %>%
    summarize(Cell_Count = min(Cell_Count), .groups = "drop"),
  
  "Z14004" = clone_data %>%
    filter(Subject == "Z14004") %>%
    group_by(Time) %>%
    summarize(Cell_Count = min(Cell_Count), .groups = "drop")
)

sample_sizes_vector <- list(
  "Z13264" = data_sample_sizes[["Z13264"]]$Cell_Count,
  "Z14004" = data_sample_sizes[["Z14004"]]$Cell_Count
)
names(sample_sizes_vector[["Z13264"]]) <- paste0("t", tps[["Z13264"]])
names(sample_sizes_vector[["Z14004"]]) <- paste0("t", tps[["Z14004"]])

# Get contributions histogram breaks
break_list <- list(
  "Z13264" = lapply(tps[["Z13264"]], function(tp) {
    curr_data <- clone_data %>%
      filter(Subject == "Z13264" & Time == tp)
    
    hist_data <- hist(log10(curr_data$Clone_Size/curr_data$Cell_Count), plot = FALSE)
    
    return(c(-Inf, hist_data$breaks, 0))
  }),
  
  "Z14004" = lapply(tps[["Z14004"]], function(tp) {
    curr_data <- clone_data %>%
      filter(Subject == "Z14004" & Time == tp)
    
    hist_data <- hist(log10(curr_data$Clone_Size/curr_data$Cell_Count), plot = FALSE)
    
    return(c(-Inf, hist_data$breaks, 0))
  })
)

names(break_list[["Z13264"]]) <- paste0("t", tps[["Z13264"]])
names(break_list[["Z14004"]]) <- paste0("t", tps[["Z14004"]])

get_contribution_dist_df_local <- function(sim_data, subject = "Z13264") {
  return(data.frame(do.call(rbind, lapply(tps[[subject]], function(tp) {
    curr_data <- sim_data[, paste0("t", tp)]
    
    existing_clones <- which(curr_data > 0)
    existing_clone_sizes <- curr_data[existing_clones]
    
    cell_count <- sum(existing_clone_sizes)
    hist_data <- hist(
      log10(existing_clone_sizes/cell_count), 
      plot = FALSE, 
      breaks = break_list[[subject]][[paste0("t", tp)]]
    )
    non_zero <- which(hist_data$counts > 0)
    
    return(data.frame(
      Time = tp,
      
      Clone_Contribution = 10**hist_data$breaks[non_zero],
      Frequency = hist_data$counts[non_zero]/cell_count
    ))
  }))))
}

# Run simulation ####
run_sim_and_sample <- function(
    nr_clones,
    pA, dA, tQA, tAQ, pQ,
    kA, kQ,
    init_A, init_Q,
    subject
) {
  sim_data <- simulate_HSC_simplified(
    nr_clones,
    pA, dA, tQA, tAQ, pQ,
    kA, kQ,
    init_A, init_Q,
    tps[[subject]]
  )
  
  colnames(sim_data$Active) <- paste0("t", c(tps[[subject]]))
  
  sampled_data <- sample_sim_data(
    sim_data$Active,
    tps[[subject]],
    sample_sizes_vector[[subject]]
  )
  
  return(sampled_data)
}


# Get differences ####
get_differences <- function(sim_results, nr_sims, subject) {
  # Diversity: ####
  diversity_df <- 
    # Get the mean simulation results
    do.call(rbind, lapply(1:nr_sims, function(sim_idx) 
      get_diversity_df(sim_results[[sim_idx]])
    )) %>%
    group_by(Time) %>%
    summarize(Nr_Clones = mean(Nr_Clones), .groups = "drop") %>%
    # Compare with Data
    merge(
      clones_over_time %>% filter(Subject == subject),
      by = c("Time"),
      all = FALSE,
      suffixes = c("_Sim", "_Data")
    ) %>%
    filter(Time > 365) %>%
    mutate(
      Res = log10(Nr_Clones_Sim) - log10(Nr_Clones_Data)
    ) %>%
    # Return points and residuals
    select(
      Time,
      Nr_Clones_Sim, Nr_Clones_Data,
      Res
    )
  
  # Clone Contribution: ####
  compared_tps <- ifelse(
    rep(subject == "Z13264", 2),
    c(397, 462),
    c(362, 495)
  )
  
  contribution_df <- 
    # Get the mean simulation results
    do.call(rbind, lapply(1:nr_sims, function(sim_idx) 
      get_contribution_dist_df_local(sim_results[[sim_idx]], subject)
    )) %>%
    group_by(Time, Clone_Contribution) %>%
    summarize(Frequency = mean(Frequency), .groups = "drop") %>%
    # Compare with Data
    merge(
      clone_contribution_dist %>% filter(Subject == subject),
      by = c("Time", "Clone_Contribution"),
      all = TRUE,
      suffixes = c("_Sim", "_Data")
    ) %>%
    filter(Time %in% compared_tps & Frequency_Data > 1e-3) %>%
    # filter(Time %in% compared_tps & Frequency_Data > 1e-) %>%
    mutate(
      # Before taking difference, change NA to 1 so when log10 its applied it becomes 0
      Res = log10(coalesce(Frequency_Sim, 1)) -  log10(coalesce(Frequency_Data, 1))
    ) %>%
    # Return points and residuals
    select(
      Time, Clone_Contribution,
      Frequency_Sim, Frequency_Data,
      Res
    )
  
  # Single Occurring Clones: ####
  so_df <- 
    # Get the mean simulation results
    do.call(rbind, lapply(1:nr_sims, function(sim_idx) 
      get_so_df(sim_results[[sim_idx]])
    )) %>%
    group_by(Time) %>%
    summarize(
      SO_Clone_Contribution = mean(SO_Clone_Contribution),
      .groups = "drop"
    ) %>%
    # Compare with Data
    merge(
      so_occurence %>% filter(Subject == subject),
      by = c("Time"),
      all = FALSE,
      suffixes = c("_Sim", "_Data")
    ) %>%
    filter(Time > 365 & Total_Cells > 1e3) %>%
    mutate(
      # Before taking the log, change 0 to 1 so when log10 its applied it becomes 0
      Res =
        log10(SO_Clone_Contribution_Sim + (SO_Clone_Contribution_Sim == 0)) -
        log10(SO_Clone_Contribution_Data + (SO_Clone_Contribution_Data == 0))
    ) %>%
    # Return points and residuals
    select(
      Time,
      SO_Clone_Contribution_Sim, SO_Clone_Contribution_Data,
      Res
    )
  
  # Return dfs ####
  return(list(
    "Diversity" = diversity_df,
    "Contribution" = contribution_df,
    "SO" = so_df
  ))
}

# Optimization Functions ####
run_two_compartments <- function(
    optimized_params,
    subject = NULL,
    nr_sims = 1,
    sim_parallel = FALSE,
    return_result = FALSE
) {
  
  if(is.null(subject)) {
    # print("No subject specified. Running both.")
    
    return(list(
      "Z13264" = run_two_compartments(
        optimized_params,
        "Z13264",
        nr_sims,
        sim_parallel,
        return_result
      ),
      "Z14004" = run_two_compartments(
        optimized_params,
        "Z14004",
        nr_sims,
        sim_parallel,
        return_result
      )
    ))
  } else {
    # print(paste("Running subject", subject))
  }
  
  # Init simulation parameters ####
  pars <- base_params
  for(name in names(optimized_params)) {
    pars[[name]] <- optimized_params[[name]]
  }
  
  pA <- pars[["pA"]]
  dA <- pars[["dA"]]
  tQA <- pars[["tQA"]]
  tAQ <- pars[["tAQ"]]
  pQ <- pars[["pQ"]]
  
  nr_clones <- ceiling(pars[["clone_mult"]]*1e4)
  if(nr_clones <= 0) {
    stop("Nr clones is 0.")
  }
  
  init_cell_dist <- sample_sizes_from_contribution(
    distribution_time[[subject]],
    nr_clones,
    sbj = subject,
    use_inferred_data = TRUE
  )
  init_cell_count <- sum(init_cell_dist)
  
  if(init_cell_count <= 0) {
    stop("Init cell count is 0.")
  }
  
  init_cell_vector <- rep(1:nr_clones, init_cell_dist)
  
  # Select cells to start active or quiescent
  q_index <- sample(
    init_cell_vector,
    floor(init_cell_count*pars[["q_proportion"]])
  )
  q_clone_sizes <- table(q_index)
  init_Q <- rep(0, nr_clones)
  init_Q[as.integer(names(q_clone_sizes))] <- q_clone_sizes
  init_A <- init_cell_dist - init_Q
  
  # Set carrying capacity based on equilibrium
  eq_A <- pars[["space_mult_A"]]*sum(init_A)
  eq_Q <- pars[["space_mult_Q"]]*sum(init_Q)
  
  if(eq_A <= 0) {
    stop("Active population on equilibrium is 0")
  }
  
  if(eq_Q <= 0) {
    stop("Quiescent population on equilibrium is 0")
  }
  
  kA <- eq_A*pA/(pA - dA)
  kQ <- eq_Q*(tAQ*pA*eq_A)/(tAQ*pA*eq_A - tQA*dA*eq_Q)
  if(kA < 0) {
    stop("Negative kA")
  }
  if(kQ < 0) {
    stop("Negative kQ")
  }
  
  loop_fun <- ifelse(sim_parallel, parallel::mclapply, lapply)
  
  # Simulate ####
  sim_results <- loop_fun(1:nr_sims, function(sim_idx) 
    return(run_sim_and_sample(
      nr_clones,
      pA, dA, tQA, tAQ, pQ,
      kA, kQ,
      init_A, init_Q,
      subject
    ))
  )
  if(return_result) {
    return(sim_results)
  }
  
  # Return compared data points ####
  return(get_differences(sim_results, nr_sims, subject))
}

run_single_compartment <- function(
    optimized_params,
    subject = NULL,
    nr_sims = 1,
    sim_parallel = FALSE,
    return_result = FALSE
) {
  
  if(is.null(subject)) {
    return(list(
      "Z13264" = run_single_compartment(
        optimized_params,
        "Z13264",
        nr_sims,
        sim_parallel,
        return_result
      ),
      "Z14004" = run_single_compartment(
        optimized_params,
        "Z14004",
        nr_sims,
        sim_parallel,
        return_result
      )
    ))
  }
  
  # Init simulation parameters ####
  pars <- base_params
  for(name in names(optimized_params)) {
    pars[[name]] <- optimized_params[[name]]
  }
  
  pA <- pars[["pA"]]
  dA <- pars[["dA"]]
  
  nr_clones <- ceiling(pars[["clone_mult"]]*1e4)
  
  init_A <- sample_sizes_from_contribution(
    distribution_time[[subject]],
    nr_clones,
    sbj = subject,
    use_inferred_data = TRUE
  )
  
  # Set carrying capacity based on equilibrium
  eq_A <- pars[["space_mult_A"]]*sum(init_A)
  
  kA <- eq_A*pA/(pA - dA)
  
  loop_fun <- ifelse(sim_parallel, parallel::mclapply, lapply)
  
  # Simulate ####
  sim_results <- loop_fun(1:nr_sims, function(sim_idx)
    run_sim_and_sample(
      nr_clones,
      pA, dA, 0, 0, 0,
      kA, 1,
      init_A, 0,
      subject
    )
  )
  if(return_result) {
    return(sim_results)
  }
  
  # Return compared data points ####
  return(get_differences(sim_results, nr_sims, subject))
}
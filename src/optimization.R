load("data/derived_radtke.Rda")
source("src/utils.R")
source("src/Simplified_HSC_Model.R")

library(dplyr)

# Base optimization params (needed for single compartment model) ####
base_params = c(
  tQA = 0,
  tAQ = 0,
  pQ = 0,
  space_mult_Q = 1
)

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
      get_contribution_dist_df(
        sim_results[[sim_idx]], 
        break_list = break_list[[subject]]
      )
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
    filter(Time %in% compared_tps & Clone_Contribution <= 1e-3) %>%
    # filter(Time %in% compared_tps) %>%
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
run_one_cell <- function(
    optimized_params,
    subject = NULL,
    nr_sims = 1,
    sim_parallel = FALSE,
    return_result = FALSE
) {
  if(is.null(subject)) {
    # print("No subject specified. Running both.")
    
    return(list(
      "Z13264" = run_one_cell(
        optimized_params,
        "Z13264",
        nr_sims,
        sim_parallel,
        return_result
      ),
      "Z14004" = run_one_cell(
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
  
  init_A <- 1
  init_Q <- 0
  
  # Set carrying capacity based on equilibrium
  eq_A <- pars[["space_mult_A"]]*nr_clones
  eq_Q <- pars[["space_mult_Q"]]*nr_clones
  
  if(eq_A <= 0) {
    stop("Active population on equilibrium is 0")
  }
  
  if(eq_Q <= 0) {
    stop("Quiescent population on equilibrium is 0")
  }
  
  kA <- eq_A*pA/(pA - dA)
  
  kQ <- ifelse(
    tAQ > 0,
    eq_Q*(tAQ*pA*eq_A)/(tAQ*pA*eq_A - tQA*dA*eq_Q),
    1
  )
  if(kA < 0) {
    stop("Negative kA")
  }
  if(kQ < 0) {
    stop("Negative kQ")
  }
  
  loop_fun <- ifelse(sim_parallel, parallel::mclapply, lapply)
  
  # Simulate ####
  sim_results <- loop_fun(1:nr_sims, function(sim_idx) 
    return(simulate_HSC_with_sampling(
      nr_clones,
      pA, dA, tQA, tAQ, pQ,
      kA, kQ,
      init_A, init_Q,
      data_tps[[subject]],
      sample_sizes_vector[[subject]]
    ))
  )
  if(return_result) {
    return(sim_results)
  }
  
  # Return compared data points ####
  return(get_differences(sim_results, nr_sims, subject))
}
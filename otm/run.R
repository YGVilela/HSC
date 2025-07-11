# Execution vars ####

## Def ####
base_path <- "data"
out_one <- file.path(base_path, "one")
out_two <- file.path(base_path, "two")
out_unified <- file.path(base_path, "unified")
generationSize <- 2
generations <- 1
nr_sims <- 1
cores <- 1

## Checks ####
if(cores > parallel::detectCores()) {
  stop(paste0(
    "Trying to use more cores (", cores,
    ") than you have... (", parallel::detectCores(),
    ")."
  ))
}

dir.create(out_one, recursive = TRUE, showWarnings = FALSE)
dir.create(out_two, recursive = TRUE, showWarnings = FALSE)
dir.create(out_unified, recursive = TRUE, showWarnings = FALSE)


# External libraries ####
list.of.packages <- c("dplyr", "GA", "parallel", "doParallel")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)
library(GA)

# Data ####
load("weights.Rda")
load("radtke_data.Rda")

# Simulation function ####
simulate_HSC <- function(
    simulation_pars,
    init_A, init_Q,
    tps,
    sample_sizes = NULL
) {
  with(simulation_pars, {
    
    base_rate <- 1
    
    curr_A <- rep(1:(nr_clones), init_A)
    curr_Q <- rep(1:nr_clones, init_Q)
    
    max_time <- max(tps)
    active_history <- matrix(
      data = 0,
      nrow = nr_clones,
      ncol = length(tps),
      dimnames = list(
        paste0("x", 1:nr_clones),
        paste0("t", tps)
      )
    )
    
    for(tp in 1:max_time) {
      # Get total population
      size_A <- length(curr_A)
      size_Q <- length(curr_Q)
      
      # Update rates
      prol_rate <- pA*(1 - size_A/kA)
      diff_rate <- dA
      deact_rate <- tAQ*(1 - size_Q/kQ)
      act_rate <- tQA*(1 - size_A/kA)
      simple_act_rate <- act_rate*(1 - pQ)
      prol_act_rate <- act_rate*pQ
      
      # Check rates
      if(diff_rate+prol_rate+deact_rate > base_rate) {
        stop(paste(
          "Total rate on active compartment is higher than base rate", 
          diff_rate+prol_rate+deact_rate
        ))
      }
      if(act_rate > base_rate) {
        stop(paste(
          "Total rate on quiescent compartment is higher than base rate",
          act_rate
        ))
      }
      
      # Sample active cells
      fate_active <- runif(size_A, min = 0, max = base_rate)
      diff_idx <- which(fate_active < diff_rate)
      prol_idx <- which(fate_active >= diff_rate & fate_active < diff_rate+prol_rate)
      deact_idx <- which(fate_active >= diff_rate+prol_rate & fate_active < diff_rate+prol_rate+deact_rate)
      idle_active_idx <- which(fate_active >= diff_rate+prol_rate+deact_rate)
      
      # Sample quiescent cells
      fate_quiescent <- runif(size_Q, min = 0, max = base_rate)
      simple_act_idx <- which(fate_quiescent < simple_act_rate)
      prol_act_idx <- which(fate_quiescent >= simple_act_rate & fate_quiescent < simple_act_rate+prol_act_rate)
      idle_quiescent_idx <- which(fate_quiescent >= simple_act_rate+prol_act_rate)
      
      # Update cell pools
      next_A <- c(
        # Active cells that did nothing
        curr_A[idle_active_idx],
        # Twice the active cells that proliferated
        curr_A[prol_idx], curr_A[prol_idx],
        # Quiescent cells that activated
        curr_Q[simple_act_idx],
        # Twice the quiescent cells that proliferated and activated
        curr_Q[prol_act_idx], curr_Q[prol_act_idx]
      )
      
      next_Q <- c(
        # Quiescent cells that did nothing
        curr_Q[idle_quiescent_idx],
        # Active cells that deactivated
        curr_A[deact_idx]
      )
      
      curr_A <- next_A
      curr_Q <- next_Q
      
      if(tp %in% tps) {
        if(length(curr_A) > 0) {
          if(is.null(sample_sizes)) {
            tmp <- table(curr_A)
            
            active_history[paste0("x", names(tmp)), paste0("t", tp)] <- tmp
          } else {
            sampled_HSC <- sample(
              1:length(curr_A),
              min(sample_sizes[paste0("t", tp)], length(curr_A))
            )
            tmp <- table(curr_A[sampled_HSC])
            
            active_history[paste0("x", names(tmp)), paste0("t", tp)] <- tmp
            curr_A <- curr_A[-sampled_HSC]
          }
          
          
        }
      }
    }
    
    return(active_history)
  })
}

# Metric evaluation ####
get_diversity_df <- function(sim_data) {
  tps <- as.integer(substring(colnames(sim_data), 2))
  
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
    Nr_Clones = nr_clones
  ))
}

get_contribution_dist_df <- function(sim_data, break_list = NULL, keep_all = TRUE) {
  tps <- as.integer(substring(colnames(sim_data), 2))
  
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
    kept_rows <- which(hist_data$breaks > -Inf & hist_data$breaks < 0)
    
    return(data.frame(
      Time = tp,
      Clone_Contribution = 10**hist_data$breaks[kept_rows],
      Frequency = hist_data$counts[kept_rows]/length(existing_clones)
    ))
  }))))
}

get_so_df <- function(sim_data) {
  tps <- as.integer(sapply(colnames(sim_data), function(t_name) substring(t_name, 2)))
  
  total_cells <- colSums(sim_data[, paste0("t", tps)])
  nr_clones <- colSums(sim_data[, paste0("t", tps)] > 0)
  
  clone_occurrences <- rowSums(sim_data[, paste0("t", tps)] > 0)
  so_clones <- which(clone_occurrences == 1)
  nr_so_clones <- colSums(sim_data[so_clones, paste0("t", tps)] > 0)
  
  return(data.frame(
    Total_Cells = total_cells,
    SO_Clone_Contribution = nr_so_clones/nr_clones,
    Time = tps
  ))
}

# Optimization routine ####
transform_pars <- function(otm_pars) {
  base_params = c(
    tQA = 0,
    tAQ = 0,
    pQ = 0,
    space_mult_Q = 1
  )
  
  aux_pars <- base_params
  for(name in names(otm_pars)) {
    aux_pars[[name]] <- otm_pars[[name]]
  }
  
  pA <- aux_pars[["pA"]]
  dA <- aux_pars[["dA"]]
  tQA <- aux_pars[["tQA"]]
  tAQ <- aux_pars[["tAQ"]]
  pQ <- aux_pars[["pQ"]]
  
  nr_clones <- ceiling(aux_pars[["clone_mult"]]*1e4)
  if(nr_clones <= 0) {
    stop("Nr clones is 0.")
  }
  
  init_A <- 1
  init_Q <- 0
  
  # Set carrying capacity based on equilibrium
  eq_A <- aux_pars[["space_mult_A"]]*nr_clones
  eq_Q <- aux_pars[["space_mult_Q"]]*nr_clones
  
  if(eq_A <= 0) {
    stop("Active population on equilibrium is 0")
  }
  
  if(eq_Q <= 0) {
    stop("Quiescent population on equilibrium is 0")
  }
  
  kA <- eq_A*pA/(pA - dA)
  
  kQ <- ifelse(
    tAQ > 0 & tQA > 0,
    eq_Q*(tAQ*pA*eq_A)/(tAQ*pA*eq_A - tQA*dA*eq_Q),
    1
  )
  if(kA < 0) {
    stop("Negative kA")
  }
  if(kQ < 0) {
    stop("Negative kQ")
  }
  
  return(list(
    nr_clones = nr_clones,
    
    pA = pA,
    dA = dA,
    tAQ = tAQ,
    tQA = tQA,
    pQ = pQ,
    
    kA = kA,
    kQ = kQ
  ))
}

evaluate_res <- function(sim_results, animal_id) {
  eps <- 1e-8
  # Diversity: ####
  diversity_df <- 
    get_diversity_df(sim_results) %>%
    # Compare with Data
    merge(
      clones_over_time %>% filter(Animal_Id == animal_id),
      by = c("Time"),
      all = FALSE,
      suffixes = c("_Sim", "_Data")
    ) %>%
    filter(Time > 365) %>%
    mutate(
      Res = log10(coalesce(Nr_Clones_Sim, 0) + eps) - log10(coalesce(Nr_Clones_Data, 0) + eps)
    ) %>%
    # Return points and residuals
    mutate(Graph = "Diversity") %>%
    select(Graph, Time, Nr_Clones_Sim, Res) %>% 
    rename(
      X = Time,
      Y = Nr_Clones_Sim
    )
  
  # Clone Contribution: ####
  compared_tps <- ifelse(
    rep(animal_id == "Z13264", 2),
    c(397, 462),
    c(362, 495)
  )
  
  contribution_df <- 
    # Get the mean simulation results
    get_contribution_dist_df(
      sim_results, 
      break_list = break_list[[animal_id]]
    )  %>%
    # Compare with Data
    merge(
      clone_contribution_dist %>% filter(Animal_Id == animal_id),
      by = c("Time", "Clone_Contribution"),
      all = TRUE,
      suffixes = c("_Sim", "_Data")
    ) %>%
    filter(
      Clone_Contribution <= 1e-3 &
      Time %in% compared_tps
    ) %>%
    mutate(
      Res = log10(coalesce(Frequency_Sim, 0) + eps) -  log10(coalesce(Frequency_Data, 0) + eps)
    ) %>%
    # Return points and residuals
    mutate(Graph = paste0("Contribution_", Time)) %>%
    select(Graph, Clone_Contribution, Frequency_Sim, Res) %>%
    rename(
      Y = Frequency_Sim,
      X = Clone_Contribution
    )
  
  # Single Occurring Clones: ####
  so_df <- 
    # Get the mean simulation results
    get_so_df(sim_results) %>%
    # Compare with Data
    merge(
      so_occurence %>% filter(Animal_Id == animal_id),
      by = c("Time"),
      all = FALSE,
      suffixes = c("_Sim", "_Data")
    ) %>%
    filter(Time > 365 & Total_Cells_Data > 1e2) %>%
    mutate(
      Res =
        log10(coalesce(SO_Clone_Contribution_Sim, 0) + eps) -
        log10(coalesce(SO_Clone_Contribution_Data, 0) +  eps)
    ) %>%
    # Return points and residuals
    mutate(Graph = "SO") %>%
    select(Graph, Time, SO_Clone_Contribution_Sim, Res) %>%
    rename(
      X = Time,
      Y = SO_Clone_Contribution_Sim
    )
  
  # Return dfs ####
  return(bind_rows(diversity_df, contribution_df, so_df))
}

evaluate_metric <- function(full_res) {
  # Merge residuals with weights df
  merge(full_res, weight_df) %>%
    # Normalize residual with graph weight
    mutate(Normalized_Res = Res*Weight_Graph) %>%
    # Evaluate mean distance (L1-norm) for each graph in each simulation
    group_by(Animal_Id, Graph, Sim_Idx) %>%
    summarise(Mean_L1 = mean(abs(Normalized_Res)), .groups = "drop") %>%
    # Evaluate mean distance for each graph across simulations
    group_by(Animal_Id, Graph) %>%
    summarise(Mean_L1 = mean(Mean_L1), .groups = "drop") %>%
    # Sum the distance of all graph
    select(Mean_L1) %>%
    sum() %>%
    return()
}

obj_func <- function(otm_pars, par_names, nr_sims) {
  names(otm_pars) <- par_names
  simulation_pars <- transform_pars(otm_pars)
  
  res_df_new <- lapply(1:nr_sims, function(sim_idx) {
    animal_id_list <- c("Z14004", "Z13264")
    
    residual_list <- lapply(animal_id_list, function(animal_id)
      simulate_HSC(
        simulation_pars,
        init_A = 1, init_Q = 0,
        tps = data_tps[[animal_id]],
        sample_sizes = sample_sizes_vector[[animal_id]]
      ) %>% 
        evaluate_res(animal_id) %>%
        mutate(Animal_Id = animal_id)
    )
    
    return(
      bind_rows(residual_list) %>%
      mutate(Sim_Idx = sim_idx)
    )
  }) %>% bind_rows()
  
  obj_value <- evaluate_metric(bind_rows(res_df_new))
  
  return(obj_value)
}

# GA ####

postfit <- function(curr_ga, base_iter_path) {
  iter <- curr_ga@iter
  
  # save info
  print(Sys.time())
  save(
    curr_ga,
    file = paste0(base_iter_path, iter, ".Rda")
  )
  
  curr_ga 
}

run_ga <- function(
    bounds,
    base_iter_path,
    final_path,
    generationSize,
    generations,
    nr_sims,
    cores = 2
  ) {
  
  start_time <- Sys.time()
  ga_res <- ga(
    type = "real-valued",
    fitness = function(x) tryCatch(
      -obj_func(x, bounds$names, nr_sims),
      error = function(e) {
        print(conditionMessage(e))
        
        return(-1e10)
      }
    ),
    names = bounds$names,
    lower = bounds$lower,
    upper = bounds$upper,
    popSize = generationSize,
    maxiter = generations,
    run = 10,
    parallel = cores,
    postFitness = function(obj, ...) postfit(obj, base_iter_path),
    monitor = TRUE,
    keepBest = TRUE
  )
  end_time <- Sys.time()
  print(end_time - start_time)
  
  save(
    ga_res,
    file = final_path
  )
  print(paste("Final result saved on", final_path))
}

# Run! ####
## One Compartment ####
bounds_one <- list(
  "names" = c("pA", "dA", "clone_mult", "space_mult_A"),
  "lower" = c(0.4,  0.025,  12,  1),
  "upper" = c(0.9, 0.15, 50, 5)
)

run_ga(
  bounds = bounds_one,
  base_iter_path = file.path(out_one, "iter_"),
  final_path = file.path(out_one, "final.Rda"),
  generationSize = generationSize,
  generations = generations,
  nr_sims = nr_sims,
  cores = cores
)

## Two Compartments ####
bounds_two <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q"),
  "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1),
  "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5)
)

run_ga(
  bounds = bounds_two,
  base_iter_path = file.path(out_two, "iter_"),
  final_path = file.path(out_two, "final.Rda"),
  generationSize = generationSize,
  generations = generations,
  nr_sims = nr_sims,
  cores = cores
)

## Unified ####
bounds_unified <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ"),
  "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1, -0.05),
  "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5, 1)
)

run_ga(
  bounds = bounds_unified,
  base_iter_path = file.path(out_unified, "iter_"),
  final_path = file.path(out_unified, "final.Rda"),
  generationSize = generationSize,
  generations = generations,
  nr_sims = nr_sims,
  cores = cores
)
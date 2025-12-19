# Execution vars ####

## Def ####
base_path <- "last_run"
out_unified <- file.path(base_path, "unified")
out_two <- file.path(base_path, "two")
generationSize <- 1000
generations <- 50
nr_sims <- 10
cores <- 100

## Checks ####
if(cores > parallel::detectCores()) {
  stop(paste0(
    "Trying to use more cores (", cores,
    ") than you have... (", parallel::detectCores(),
    ")."
  ))
}

dir.create(out_unified, recursive = TRUE, showWarnings = FALSE)
dir.create(out_two, recursive = TRUE, showWarnings = FALSE)


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
      prol_rate <- max(pA*(1 - size_A/kA), 0)
      diff_rate <- max(dA, 0)
      deact_rate <- max(tAQ*(1 - size_Q/kQ), 0)
      act_rate <- max(tQA*(1 - size_A/kA), 0)
      simple_act_rate <- max(act_rate*(1 - pQ), 0)
      prol_act_rate <- max(act_rate*pQ, 0)
      
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
  
  pA <- max(aux_pars[["pA"]], 0)
  dA <- max(aux_pars[["dA"]], 0)
  tQA <- max(aux_pars[["tQA"]], 0)
  tAQ <- max(aux_pars[["tAQ"]], 0)
  pQ <- max(aux_pars[["pQ"]], 0)
  
  nr_clones <- ceiling(aux_pars[["clone_mult"]]*1e4)
  if(nr_clones <= 0) {
    stop("Nr clones is 0.")
  }
  
  init_A <- 1
  init_Q <- 0
  
  # Set FIXED carrying capacities
  if(tAQ > 0 & tQA > 0) {
    kA <- 1e5
    kQ <- 1e5
  } else {
    kA <- 2e5
    kQ <- 1
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
    filter(Time > 360) %>%
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
      Clone_Contribution <= 1e-3 & Time %in% compared_tps
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
    filter(Time > 360 & Total_Cells_Data > 1e2) %>%
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
    mutate(Normalized_Res = Res*Graph_Weight) %>%
    # Evaluate mean distance (Squared distance) for each graph in each simulation
    group_by(Animal_Id, Graph, Sim_Idx) %>%
    summarise(RSM = mean(Normalized_Res**2), .groups = "drop") %>%
    # Evaluate mean distance for each graph across simulations
    group_by(Animal_Id, Graph) %>%
    summarise(RSM = mean(RSM), .groups = "drop") %>%
    # Sum the distance of all graph
    select(RSM) %>%
    sum() %>%
    return()
}

obj_func <- function(otm_pars, par_names, nr_sims) {
  names(otm_pars) <- par_names
  simulation_pars <- transform_pars(otm_pars)
  
  res_df <- lapply(1:nr_sims, function(sim_idx) {
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
      bind_rows(residual_list) %>% mutate(Sim_Idx = sim_idx)
    )
  }) %>% bind_rows()
  
  obj_value <- evaluate_metric(res_df)
  
  return(obj_value)
}

# GA ####

# Save iteration data
postfit <- function(curr_ga, base_iter_path) {
  iter <- curr_ga@iter
  
  # save info
  print(Sys.time())
  save(
    curr_ga,
    file = paste0(base_iter_path, iter, ".Rda")
  )
  
  bestF <- -max(curr_ga@fitness)
  medianF <- sort(-curr_ga@fitness)[floor(curr_ga@popSize/4)]
  mutationFactor <- bestF/medianF
  save(
    mutationFactor,
    file = paste0(base_iter_path, iter, "_aux.Rda")
  )
  
  curr_ga 
}

# Mutate pQ preferentially
mutationFunc <- function(object, parent) {
  sampleProb <- rep(1, 5)
  if(length(object@names) == 6) {
    # pQ mutates twice as much
    sampleProb <- c(sampleProb, 2)
  }
  sampleProb <- sampleProb/sum(sampleProb)
  
  mutate <- as.vector(object@population[parent,])
  
  j <- sample(1:length(sampleProb), size = 1, prob = sampleProb)
  
  # When closer to fit, mutate less
  iter <- object@iter
  load(paste0(base_iter_path, iter, "_aux.Rda"))
  if(mutationFactor >= 0.5) {
    # Mutation restricted to 1% of current value
    mutate[j] <- mutate[j]+runif(1, -0.01, 0.01)*(object@upper[j] - object@lower[j])
    if(mutate[j] < object@lower[j]) {
      mutate[j] <- object@lower[j]
    }
    if(mutate[j] > object@upper[j]) {
      mutate[j] <- object@upper[j]
    }
  } else {
    mutate[j] <- runif(1, object@lower[j], object@upper[j])
  }
  
  return(mutate)
}

# When closer to fit, mutate less, mutate more often
pmutationFunc <- function(object) {
  iter <- object@iter
  load(paste0(base_iter_path, iter, "_aux.Rda"))
  if(mutationFactor >= 0.5) {
    return(0.5)
  } else {
    return(0.1)
  }
}

run_ga <- function(
    bounds,
    generationSize,
    generations,
    nr_sims,
    cores = 2,
    initial_population = NULL
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
    run = 15,
    parallel = cores,
    postFitness = function(obj, ...) postfit(obj, base_iter_path),
    monitor = TRUE,
    keepBest = TRUE,
    suggestions = initial_population,
    mutation = mutationFunc, # pQ mutates more often, mutation range decreases when closer to fixation
    pmutation = pmutationFunc, # Increase mutation rate when closer to fixation
    elitism = 1 # Only top individual guaranteed to survive
  )
  end_time <- Sys.time()
  print(end_time - start_time)
  
  save(
    ga_res,
    file = final_path
  )
  print(paste("Final result saved on", final_path))
}

# Starting with the best parameters for the two compartment as suggestions changing pQ
load("unified_suggestions.Rda")
nrSuggestions <- 100
nrMutations <- 10

ord <- order(ga_res@fitness, decreasing = TRUE)
PopSorted <- ga_res@population[ord,,drop=FALSE]
ga_res@population <- PopSorted
u <- which(!duplicated(PopSorted, margin = 1))
topInd <- PopSorted[u[1:nrSuggestions],]

## Two compartment ####
mutateInit <- function(originalPop) {
  n <- length(originalPop[,1])
  allJ <- sample(1:5, size = n, replace = TRUE)
  
  for(row in 1:n) {
    j <- allJ[row]
    newValue <- originalPop[row, j] + runif(1, -0.01, 0.01)*(ga_res@upper[j] - ga_res@lower[j])
    if(newValue < ga_res@lower[j]) {
      newValue <- ga_res@lower[j]
    }
    if(newValue > ga_res@upper[j]) {
      newValue <- ga_res@upper[j]
    }
    
    originalPop[row, j] <- newValue
  }
  
  return(originalPop)
}
initial_population <- do.call("rbind", 
  lapply(1:nrMutations, function(idx_m) mutateInit(topInd))
)
colnames(initial_population) <- c("pA", "dA", "tQA", "tAQ", "clone_mult")

bounds_two <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult"),
  "lower" = c(0.4,  0.025,  0.0001,  0.01,  3),
  "upper" = c(0.9, 0.15, 0.03, 0.075, 12)
)

base_iter_path <- file.path(out_two, "iter_")
final_path <- file.path(out_two, "final.Rda")
run_ga(
  bounds = bounds_two,
  generationSize = generationSize,
  generations = generations,
  nr_sims = nr_sims,
  cores = cores,
  initial_population = initial_population
)

## Unified #### 
pQRange <- 0.1*(0:(nrMutations - 1))/(nrMutations-1)

initial_population <- do.call("rbind", 
  lapply(pQRange, function(pQ) cbind(topInd, rep(pQ, nrSuggestions)))
)
colnames(initial_population) <- c("pA", "dA", "tQA", "tAQ", "clone_mult", "pQ")

bounds_unified <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "pQ"),
  "lower" = c(0.4,  0.025,  0.0001,  0.01,  3, -0.001),
  "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 0.5)
)

base_iter_path <- file.path(out_unified, "iter_")
final_path <- file.path(out_unified, "final.Rda")
run_ga(
  bounds = bounds_unified,
  generationSize = generationSize,
  generations = generations,
  nr_sims = nr_sims,
  cores = cores,
  initial_population = initial_population
)

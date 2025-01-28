load("data/derived_radtke.Rda")

# Cell action sampling functions ####
sample_multinomial <- function(size, rates) {
  total_rate <- sum(rates)
  if(total_rate > 1) {
    stop(
      paste0(
        "Total rate on sample_multinomial must be smaller at most 1 (evaluated as ",
        total_rate,
        "). Try using a higher base_rate.")
    )
  }
  
  return(stats::rmultinom(
    n = 1,
    size = size,
    prob = append(rates, 1 - total_rate)
  )[1:length(rates), 1])
}

sample_deterministic <- function(size, rates) {
  return(size*rates)
}

#' @exportMethod
sample_cells <- function(
    nr_clones,
    population,
    actions,
    rates,
    sampling_method = "sample_multinomial"
) {
  # Group actions by source_compartment
  actions_by_source = list()
  for(action in actions) {
    source_compartment <- action$source_compartment
    
    actions_by_source[[source_compartment]] <- 
      append(actions_by_source[[source_compartment]], action$name)
  }
  
  # For each source compartment, use rates to pool how many cells perform
  # each action.
  sample_results <- array(
    data = 0,
    dim = c(nr_clones, length(actions)),
    dimnames = list(NULL, names(actions))
  )
  for(source_compartment in names(actions_by_source)) {
    action_list <- actions_by_source[[source_compartment]]
    # Keep for safety, but this if shouldn't be reached...
    if(length(action_list) == 0) {
      next
    }
    
    # TODO: Consider parallelizing this if there are too many clones...
    for(clone_idx in 1:nr_clones) {
      if(population[[clone_idx, source_compartment]] == 0) {
        next
      } 
      
      curr_rates <- rates[clone_idx, action_list]
      sample_results[clone_idx, action_list] <- do.call(
        sampling_method,
        list(
          size = population[clone_idx, source_compartment],
          rates = curr_rates
        )
      )
    }
  }
  
  return(sample_results)
}

# Initial size sampling functions ####
#' @exportMethod
get_sizes_from_inferred_data <- function(tp, sbj = "Z13264") {
  data <- inferred_clone_data[[sbj]][, paste0("t", tp)]
  data <- data[which(data > 0)]
  names(data) <- c()
  
  return(data)
}

#' @exportMethod
sample_sizes_from_contribution <- function(tp, nr_clones, sbj = "Z13264", use_inferred_data = FALSE) {
  if(use_inferred_data) {
    data <- inferred_clone_data[[sbj]][, paste0("t", tp)]
    data <- data[which(data > 0)]
    names(data) <- c()
  } else {
    data <- (clone_data %>% filter(Subject == sbj & Time == tp))$Clone_Size
  }
  
  contribution <- data/sum(data)
  den <- density(contribution)
  sampled_contribution <- sample(x = den$x, nr_clones, prob = den$y, replace=TRUE) + rnorm(nr_clones, 0, den$bw)
  sampled_contribution <- sampled_contribution/sum(sampled_contribution)
  sampled_sizes <- ceiling(sampled_contribution*nr_clones)
  sampled_sizes[which(sampled_sizes <= 0)] <- 1
  
  return(sampled_sizes)
}

# Simulation results sampling function ####
#' @exportMethod
sample_sim_data <- function(sim_data, tps, sample_size) {
  nr_clones <- length(sim_data[, 1])
  clone_names <- paste0("x", 1:nr_clones)
  
  sampled_data <- matrix(
    data = 0,
    nrow = nr_clones,
    ncol = length(tps),
    dimnames = list(
      clone_names,
      paste0("t", tps)
    )
  )
  
  if(length(sample_size) == 1) {
    sample_size <- rep(sample_size, length(tps))
    names(sample_size) <- paste0("t", tps)
  }
  
  for(tp in paste0("t", tps)) {
    curr_sample <- table(sample(
      rep(clone_names, sim_data[, tp]),
      min(sample_size[tp], sum(sim_data[, tp])),
      replace = FALSE
    ))
    
    sampled_data[names(curr_sample), tp] <- curr_sample
  }
  
  return(sampled_data)
}
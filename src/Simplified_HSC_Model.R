simulate_HSC <- function(
    nr_clones,
    pA, dA, tQA, tAQ, pQ,
    kA, kQ,
    init_A, init_Q,
    tps,
    sample_sizes = NULL
) {
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
  
  pb <- progress::progress_bar$new(
    format = "[:bar :current/:total] (:eta)",
    total = max_time
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
    
    pb$tick()
  }
  
  pb$terminate()
  
  return(active_history)
}

get_action_history <- function(
    nr_clones,
    pA, dA, tQA, tAQ, pQ,
    kA, kQ,
    init_A, init_Q,
    max_time
) {
  base_rate <- 1
  
  curr_A <- rep(1:(nr_clones), init_A)
  curr_Q <- rep(1:nr_clones, init_Q)
  
  action_history <- matrix(
    data = 0,
    nrow = 7,
    ncol = max_time,
    dimnames = list(
      c(
        "A", "Q",
        "Diff", "Prol", "Deact",
        "Act_Simple", "Act_Prol"
      ),
      paste0("t", 1:max_time)
    )
  )
  
  pb <- progress::progress_bar$new(
    format = "[:bar :current/:total] (:eta)",
    total = max_time
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
    
    action_history[, paste0("t", tp)] <- c(
      "A" = size_A,
      "Q" = size_Q,
      "Diff" = length(diff_idx),
      "Prol" = length(prol_idx),
      "Deact" = length(deact_idx),
      "Act_Simple" = length(simple_act_idx),
      "Act_Prol" = length(prol_act_idx)
    )
      
    
    pb$tick()
  }
  
  pb$terminate()
  
  return(action_history)
}
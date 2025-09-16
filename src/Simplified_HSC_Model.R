simulate_HSC <- function(
    nr_clones,
    pA, dA, tQA, tAQ, pQ,
    kA, kQ,
    init_A, init_Q,
    tps,
    sample_sizes = NULL,
    return_quiescent = FALSE,
    base_rate = 1
) {
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
  if(return_quiescent) {
    quiescent_history <- matrix(
      data = 0,
      nrow = nr_clones,
      ncol = length(tps),
      dimnames = list(
        paste0("x", 1:nr_clones),
        paste0("t", tps)
      )
    )
  }
  
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
      if(return_quiescent & length(curr_Q > 0)) {
        tmp <- table(curr_Q)
        
        quiescent_history[paste0("x", names(tmp)), paste0("t", tp)] <- tmp
      }
    }
    
    pb$tick()
  }
  
  pb$terminate()
  
  if(return_quiescent) {
    return(matrix(
      c(
        tps,
        t(quiescent_history),
        t(active_history)
      ),
      byrow = FALSE,
      ncol = 1+2*nr_clones,
      dimnames = list(
        NULL,
        c(
          "step",
          paste0("Q", 1:nr_clones),
          paste0("A", 1:nr_clones)
        )
      )
    ))
  } else {
    return(active_history)
  }
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

derive_analitic_times <- function(rates, kA, kQ) {
  return(with(as.list(rates), {
    
    # eq_A <- init_A
    # eq_Q <- init_Q
    # 
    # kA <- eq_A*pA/(pA - dA)
    # 
    # kQ <- eq_Q*(tAQ*pA*eq_A)/(tAQ*pA*eq_A - tQA*dA*eq_Q)
    
    mean_diff <- 1/dA
    mean_prol <- 1/dA
    mean_deact <- 1/tAQ + (kA/kQ)*(pA/(dA*tQA) - 1/tQA)
    
    mean_act <- pA/(dA*tQA)
    
    return(c(
      "Act" = mean_act,
      "Deact" = mean_deact,
      "Diff" = mean_diff,
      "Prol" = mean_prol
      # "Total Active" = 1/(1/mean_diff + 1/mean_prol + 1/mean_deact)
    ))
  }))
}

measure_empiric_times <- function(max_time, kA, kQ, rates, base_rate = 1, time_points = NULL) {
  if(is.null(time_points)) {
    time_points = c(max_time)
  }
  
  return(with(as.list(rates), {
    init_A <- kA*(pA - dA)/pA
    init_Q <- kQ*(tAQ*(pA - dA)*kA)/(tAQ*(pA - dA)*kA + tQA*dA*kQ)
    
    options("scipen"=100, "digits"=4)
    
    measured_times <- matrix(
      data = 0,
      nrow = 4,
      ncol= length(time_points),
      dimnames = list(
        c("Act", "Deact", "Diff", "Prol"),
        paste0("t", time_points)
      )
    )
    
    curr_A <- matrix(
      data = 0,
      nrow = 4,
      ncol= init_A,
      dimnames = list(
        c("Act", "Deact", "Diff", "Prol"),
        NULL
      )
    )
    curr_Q <- matrix(
      data = 0,
      nrow = 4,
      ncol= init_Q,
      dimnames = list(
        c("Act", "Deact", "Diff", "Prol"),
        NULL
      )
    )
    action_counter <- matrix(
      data = 0,
      nrow = 4,
      ncol= max_time,
      dimnames = list(
        c("Act", "Deact", "Diff", "Prol"),
        paste0("t", 1:max_time)
      )
    )
    
    # eq_A <- init_A
    # eq_Q <- init_Q
    
    # kA <- eq_A*pA/(pA - dA)
    # 
    # kQ <- eq_Q*(tAQ*pA*eq_A)/(tAQ*pA*eq_A - tQA*dA*eq_Q)
    
    pb <- progress::progress_bar$new(
      format = "[:bar :current/:total] (:eta)",
      total = max_time
    )
    for(tp in 1:max_time) {
      # Get total population
      size_A <- length(curr_A["Prol", ])
      size_Q <- length(curr_Q["Prol", ])
      
      # Update rates
      prol_rate <- pA*(1 - size_A/kA)
      diff_rate <- dA
      deact_rate <- tAQ*(1 - size_Q/kQ)
      act_rate <- tQA*(1 - size_A/kA)
      
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
      fate_active <- rmultinom(1, size_A, c(
        "diff" = diff_rate,
        "prol" = prol_rate,
        "deact" = deact_rate,
        "idle" = base_rate - (diff_rate + prol_rate + deact_rate)
      ))[,1]
      
      available_idx  <- 1:size_A
      diff_idx<- sample(available_idx , fate_active ["diff"])
      
      available_idx <- setdiff(available_idx, diff_idx)
      prol_idx <- sample(available_idx , fate_active ["prol"])
      
      available_idx <- setdiff(available_idx, prol_idx)
      deact_idx <- sample(available_idx , fate_active ["deact"])
      
      idle_active_idx <- setdiff(available_idx, deact_idx)
      
      # Sample quiescent cells
      fate_quiescent <- rbinom(1, size_Q, act_rate/base_rate)
      
      available_idx  <- 1:size_Q
      act_idx <- sample(available_idx, fate_quiescent)
      idle_quiescent_idx <- setdiff(available_idx, act_idx)
      
      ######### Compute Times
      prol_daughter <- c()
      # Prol time
      if(length(prol_idx) > 0) {
        prol_times <- tp - curr_A["Prol", prol_idx]
        temp <- table(prol_times)
        action_counter["Prol", paste0("t", names(temp))] <- action_counter["Prol", paste0("t", names(temp))] + temp
        
        prol_daughter <- curr_A[, prol_idx]
        prol_daughter["Diff", ] <- tp
        
        curr_A["Prol", prol_idx] <- tp
      }
      
      # Diff time
      if(length(diff_idx) > 0) {
        diff_times <- tp - curr_A["Diff", diff_idx]
        temp <- table(diff_times)
        action_counter["Diff", paste0("t", names(temp))] <- action_counter["Diff", paste0("t", names(temp))] + temp
        
        # No need to update timers, as these cells will be removed from the simulation
        # curr_A["Diff", diff_idx] <- tp 
      }
      
      
      # Deact time
      if(length(deact_idx) > 0) {
        deact_times <- tp - curr_A["Deact", deact_idx]
        temp <- table(deact_times)
        action_counter["Deact", paste0("t", names(temp))] <- action_counter["Deact", paste0("t", names(temp))] + temp
        
        curr_A["Deact", deact_idx] <- tp
        
        # Update quiescent timer
        curr_A["Act", deact_idx] <- curr_A["Act", deact_idx] + (tp - curr_A["Act", deact_idx])
      }
      
      # Act time
      if(length(act_idx) > 0) {
        act_times <- tp - curr_Q["Act", act_idx]
        temp <- table(act_times)
        action_counter["Act", paste0("t", names(temp))] <- action_counter["Act", paste0("t", names(temp))] + temp
        
        curr_Q["Act", act_idx] <- tp
        
        # Update active timers
        curr_Q["Prol", act_idx] <- curr_Q["Prol", act_idx] + (tp - curr_Q["Deact", act_idx])
        curr_Q["Diff", act_idx] <- curr_Q["Diff", act_idx] + (tp - curr_Q["Deact", act_idx])
        curr_Q["Deact", act_idx] <- curr_Q["Deact", act_idx] + (tp - curr_Q["Deact", act_idx])
      }
      
      ######### End Compute Times
      
      # Update cell pools
      next_A <- matrix(
        c(
          # Keep timer:
          # Active cells that did nothing
          curr_A[, idle_active_idx],
          
          # Proliferation: reset prol timer on one, the other is new
          curr_A[, prol_idx], 
          prol_daughter,
          
          # Cells that reactivated
          curr_Q[, act_idx]
        ),
        nrow = 4,
        dimnames = list(
          c("Act", "Deact", "Diff", "Prol"),
          NULL
        ),
        byrow = FALSE
      )
      
      next_Q <- matrix(
        c(
          # Keep timer:
          # Quiescent cells that did nothing
          curr_Q[, idle_quiescent_idx],
          
          # Cells that deactivated
          curr_A[, deact_idx]
        ),
        nrow = 4,
        dimnames = list(
          c("Act", "Deact", "Diff", "Prol"),
          NULL
        ),
        byrow = FALSE
      )
      
      curr_A <- next_A
      curr_Q <- next_Q
      
      if(tp %in% time_points) {
        measured_times[, paste0("t", tp)] <- rowSums(t(
          t(action_counter[, 1:tp])*1:tp
        ))/rowSums(action_counter[, 1:tp])
      }
      
      pb$tick()
    }
    
    pb$terminate()
    
    options("scipen"=0, "digits"=7)
    return(list(
      Counter = action_counter,
      Measured_Times = measured_times
    ))
  }))
}

solve_with_ODE <- function(
    nr_clones,
    pA, dA, tQA, tAQ, pQ,
    kA, kQ,
    init_A, init_Q,
    tps
) {
  vector_field <- function(t, state, ...) {
    quiesc_sum <- sum(sapply(
      paste0("Q", 1:nr_clones), 
      function(clone_id) state[[clone_id]]
    ))
    active_sum <- sum(sapply(
      paste0("A", 1:nr_clones), 
      function(clone_id) state[[clone_id]]
    ))
    
    quiesc_field <- sapply(1:nr_clones, function(clone_idx) {
      Q <- state[[paste0("Q", clone_idx)]]
      A <- state[[paste0("A", clone_idx)]]
      
      return(
        -tQA*Q*(1 - active_sum/kA) + tAQ*A*(1 - quiesc_sum/kQ)
      )
    })
    active_field <- sapply(1:nr_clones, function(clone_idx) {
      Q <- state[[paste0("Q", clone_idx)]]
      A <- state[[paste0("A", clone_idx)]]
      
      return(
        tQA*Q*(1 - active_sum/kA) - tAQ*A*(1 - quiesc_sum/kQ) +
          (tQA*pQ*Q + pA*A)*(1 - active_sum/kA) - dA*A
      )
    })
    
    return(
      list(c(quiesc_field, active_field))
    )
  }
  
  if(length(init_A) == 1) {
    init_A <- rep(init_A, nr_clones)
  }
  if(length(init_Q) == 1) {
    init_Q <- rep(init_Q, nr_clones)
  }
  
  x0 <- c(
    init_Q,
    init_A
  )
  names(x0) <- c(
    paste0("Q", 1:nr_clones),
    paste0("A", 1:nr_clones)
  )
  
  return(
    deSolve::ode(
      y = x0, times = tps, func = vector_field
    )
  )
}

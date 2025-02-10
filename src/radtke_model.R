#' Slight modifications from https://github.com/KiemLab-RIS/HSC-Clone-Differentiation-Mathematical-Model/blob/main/HSC%20Clone%20Differentiation%20Mathematical%20Model.R
#' @exportMethod 
simulate_radtke <- function(
    nr_clones,
    division_rate,
    tps,
    init = 1
){
  curr_cells = rep(1:nr_clones, init)
  sim_hist = array(
    data = 0,
    dim = c(
      nr_clones,
      length(tps)
    ),
    dimnames = list(
      paste0("x", 1:nr_clones),
      paste0("t", tps)
    )
  )
  
  last_day <- max(tps)
  pb <- progress::progress_bar$new(
    format = "[:bar :current/:total] (:eta)",
    total = last_day
  )
  for (day in 1:last_day) {
    # pick cells to divide
    div.index = sample(
      1:length(curr_cells),
      floor(length(curr_cells) * division_rate),
      replace = FALSE
    )
    div.cells = curr_cells[div.index]
    
    # remove div cells from current cells
    curr_cells = curr_cells[-div.index]
    
    # duplicate div cells
    ly = c(div.cells,div.cells)
    
    # fate, stay stem or differentiate
    # In order to have homeostasis, half the cells must differentiate and half
    # must proliferate.
    fate = runif(length(ly))
    new.stem.index = which(fate >= 0.5)
    new.stem = ly[new.stem.index]
    
    # add new stem cells
    curr_cells = c(curr_cells,new.stem)
    
    # Save clone data
    if(day %in% tps) {
      tmp = rle(sort(curr_cells))
      sim_hist[paste0("x", tmp$values), paste0("t", day)] = tmp$lengths
    }
    pb$tick()
  }
  
  pb$terminate()
  
  return(sim_hist)
}

simulate_radtke_with_samples <- function(
    nr_clones,
    division_rate,
    tps,
    sample_sizes,
    init = 1
){
  curr_cells = rep(1:nr_clones, init)
  sim_hist = array(
    data = 0,
    dim = c(
      nr_clones,
      length(tps)
    ),
    dimnames = list(
      paste0("x", 1:nr_clones),
      paste0("t", tps)
    )
  )
  
  last_day <- max(tps)
  pb <- progress::progress_bar$new(
    format = "[:bar :current/:total] (:eta)",
    total = last_day
  )
  for (day in 1:last_day) {
    # pick cells to divide
    div.index = sample(
      1:length(curr_cells),
      floor(length(curr_cells) * division_rate),
      replace = FALSE
    )
    div.cells = curr_cells[div.index]
    
    # remove div cells from current cells
    curr_cells = curr_cells[-div.index]
    
    # duplicate div cells
    ly = c(div.cells,div.cells)
    
    # fate, stay stem or differentiate
    # In order to have homeostasis, half the cells must differentiate and half
    # must proliferate.
    fate = runif(length(ly))
    new.stem.index = which(fate >= 0.5)
    new.stem = ly[new.stem.index]
    
    # add new stem cells
    curr_cells = c(curr_cells,new.stem)
    
    # Save clone data
    if(day %in% tps) {
      sampled_cells = sample(
        1:length(curr_cells),
        min(length(curr_cells), sample_sizes[paste0("t", day)])
      )
      
      if(length(sampled_cells > 0)) {
        tmp = table(curr_cells[sampled_cells])
      
        sim_hist[paste0("x", names(tmp)), paste0("t", day)] = tmp
        
        curr_cells = curr_cells[-sampled_cells]
      }
    }
    pb$tick()
  }
  
  pb$terminate()
  
  return(sim_hist)
}

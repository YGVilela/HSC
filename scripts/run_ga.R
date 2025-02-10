library(GA)

generationSize <- 100
generations <- 30
cores <- 10

# Functions ####
obj_func <- function(
    optimized_params,
    param_names
) {
  source("src/optimization.R", local = TRUE)
  load("data/one_cell/weights.Rda")
  
  names(optimized_params) <- param_names
  
  tryCatch( {
    sim_res <- run_one_cell(
      optimized_params = optimized_params,
      nr_sims = 10,
      sim_parallel = TRUE
    )
    
    diffs_df <- get_full_diffs_df(sim_res)
    
    normalized_res <-
      merge(
        diffs_df, point_weights,
        all.x = TRUE,
        all.y = FALSE
      ) %>%
      mutate(
        Normalized_Res = Res*coalesce(Weight, 1)
      )
    
    rss_df <- normalized_res %>%
      group_by(Subject, Graph) %>%
      summarise(RSS = sum(Normalized_Res**2), .groups = "drop") %>%
      merge(
        graph_weights,
        all.x = TRUE,
        all.y = FALSE
      ) %>%
      mutate(
        Normalized_RSS = RSS*coalesce(Weight, 1)
      )
    
    return(sum(rss_df$Normalized_RSS))
  },
  error = function(e) {
    print(conditionMessage(e))
    print(optimized_params)
    
    return(100000)
  }
  )
}

postfit <- function(curr_ga, base_path) {
  iter <- curr_ga@iter
  
  # save info
  print(Sys.time())
  save(
    curr_ga,
    file = paste0(base_path, iter, ".Rda")
  )
  
  curr_ga 
}

run_ga <- function(pars, ga_iter, ga_final) {
  # load(ga_final)
  # save(
  #   ga_res,
  #   file = paste0(ga_final, "_old_v3")
  # )
  
  # suggestions <- NULL
  # if(generationSize > ga_res@popSize) {
  #   stop("Are you increasing the population??")
  # } else {
  #   best_idx <- tail(sort(ga_res@fitness, index.return = TRUE)$ix, n = generationSize)
  #   suggestions <- ga_res@population[best_idx, ]
  # }
  
  start_time <- Sys.time()
  ga_res <- ga(
    type = "real-valued",
    fitness = function(x) -obj_func(x, pars$names),
    names = pars$names,
    lower = pars$lower,
    upper = pars$upper,
    popSize = generationSize,
    maxiter = generations,
    run = 10,
    parallel = cores,
    postFitness = function(obj, ...) postfit(obj, ga_iter),
    monitor = TRUE,
    # suggestions = suggestions,
    keepBest = TRUE
  )
  end_time <- Sys.time()
  print(end_time - start_time)
  
  save(
    ga_res,
    file = ga_final
  )
  print(paste("Final result saved on", ga_final))
}

# Two compartments ####
pars_two_comp <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q"),
  "lower" = c(0.6,  0.1,  0.001,  0.04,  4,  1,  1),
  "upper" = c(0.75, 0.15, 0.006, 0.07, 6, 3, 2.5)
  # "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1),
  # "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5)
)
ga_iter_two_comp <- "data/one_cell_v2/ga_two_"
ga_final_two_comp <- "data/one_cell_v2/ga_two_final.Rda"

run_ga(
  pars_two_comp,
  ga_iter_two_comp,
  ga_final_two_comp
)

# With pQ ####
pars_pQ <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ"),
  "lower" = c(0.6,  0.1,  0.001,  0.04,  4,  1,  1, -0.05),
  "upper" = c(0.75, 0.15, 0.006, 0.07, 6, 3, 2.5, 0.5)
  # "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1, -0.05),
  # "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5, 0.45)
)
ga_iter_pQ <- "data/one_cell_v2/ga_pQ_"
ga_final_pQ <- "data/one_cell_v2/ga_pQ_final.Rda"

run_ga(
  pars_pQ,
  ga_iter_pQ,
  ga_final_pQ
)

# One compartments ####
# pars_one_comp <- list(
#   "names" = c("pA", "dA", "clone_mult", "space_mult_A"),
#   "lower" = c(0.4,  0.025,  15,  1),
#   "upper" = c(0.9, 0.15, 50, 5)
# )
# ga_iter_one_comp <- "data/one_cell/ga_one_"
# ga_final_one_comp <- "data/one_cell/ga_one_final.Rda"
# 
# run_ga(
#   pars_one_comp,
#   ga_iter_one_comp,
#   ga_final_one_comp
# )

print("Done")
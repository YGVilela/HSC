library(GA)

generationSize <- 1000
generations <- 10
cores <- 100

# Functions ####
obj_func <- function(
    optimized_params,
    param_names
) {
  source("src/optimization.R", local = TRUE)
  load("data/one_cell/weights.Rda")
  
  names(optimized_params) <- param_names
  
  tryCatch( {
    sim_res <- run_no_pA(
      optimized_params = optimized_params,
      nr_sims = 1,
      sim_parallel = FALSE
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
  start_time <- Sys.time()
  ga_res <- ga(
    type = "real-valued",
    fitness = function(x) -obj_func(x, pars$names),
    names = pars$names,
    lower = pars$lower,
    upper = pars$upper,
    popSize = generationSize,
    maxiter = generations,
    # run = 10,
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

par_list <- list(
  "names" = c("dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ"),
  "lower" = c(0.01,  0.001,  0.3,  3,  0.5,  0.5, 0),
  "upper" = c( 0.1,   1,  0.9, 10,    5,    5, 1)
  # "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1),
  # "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5)
)
iter_file <- "data/no_pA/ga_"
final_file <- "data/no_pA/ga_final.Rda"

run_ga(
  par_list,
  iter_file,
  final_file
)

print("Done")
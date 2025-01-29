library(GA)

obj_func <- function(
    optimized_params
) {
  source("src/optimization.R", local = TRUE)
  load("data/two_comp_weights.Rda")
  
  names(optimized_params) <- c(
    "pA",
    "dA",
    "tQA",
    "tAQ",
    "clone_mult",
    "space_mult_A",
    "space_mult_Q"
  )
  
  tryCatch( {
    sim_res <- run_two_compartments(
      optimized_params = optimized_params
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

postfit <- function(object, ...)
{
  iter <- object@iter
  pop <- object@population
  fit <- object@fitness
  # save info
  print(Sys.time())
  save(
    pop, fit,
    file = paste0("data/ga_two_iter_", iter, ".Rda")
  )
  # output the input ga object (this is needed!!)
  object 
}

ga_out <- "ga_res.Rda"
init_population <- NULL
if(file.exists(ga_out)) {
  load(ga_out)
  
  init_population <- ga_res@population
  rm(ga_res)
}

start_time <- Sys.time()
ga_res <- ga(
  type = "real-valued",
  fitness = function(x) -obj_func(x),
  names = c(
    "pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q"
  ),
  lower = c(
    0.5, 0.05, 0.001, 0.01, 3, 1, 1
  ),
  upper = c(
    0.9, 0.1, 0.03, 0.075, 10, 5, 5
  ),
  popSize = 1000,
  maxiter = 10,
  parallel = 50,
  postFitness = postfit,
  monitor = TRUE,
  suggestions = init_population
)
end_time <- Sys.time()
print(end_time - start_time)

save(
  ga_res,
  file = ga_out
)
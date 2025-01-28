library(rgenoud)

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


genoud_out <- "genoud_two_compartments.pro"
nr_cores <- 50

cluster <- parallel::makeCluster(nr_cores)

start_time <- Sys.time()
genoud_final <- genoud(
  fn = obj_func,
  nvars = 7,
  cluster = cluster,
  pop.size = 1000,
  max.generations = 30,
  boundary.enforcement = 2,
  project.path = genoud_out,
  Domains = matrix(
    data = c(
      # pA
      0.5, 0.9,
      # dA
      0.05, 0.1,
      # tQA
      0.001, 0.03,
      # tAQ
      0.01, 0.075,
      # clone_mult
      3, 10,
      # space_mult_A
      1, 5,
      # space_mult_Q
      1, 5
    ),
    ncol = 2,
    byrow = TRUE
  )
)

end_time <- Sys.time()
print(end_time - start_time)
parallel::stopCluster(cl = cluster)
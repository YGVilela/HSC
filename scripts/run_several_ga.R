library(GA)

# Functions ####
obj_func <- function(
    optimized_params,
    param_names,
    nr_sims = 1
) {
  source("src/optimization.R", local = TRUE)
  load("data/one_cell/weights.Rda")
  
  names(optimized_params) <- param_names
  
  tryCatch( {
    sim_res <- run_one_cell(
      optimized_params = optimized_params,
      nr_sims = nr_sims,
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
  # Initial run
  generationSize <- 1000
  generations <- 50
  cores <- 100
  start_time <- Sys.time()
  ga_res <- ga(
    type = "real-valued",
    fitness = function(x) -obj_func(x, pars$names, 1),
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
  print(paste("Broad result saved on ", ga_final))
  
  # Narrow run
  generationSize <- 100
  generations <- 10
  cores <- 100
  
  best_idx <- tail(sort(ga_res@fitness, index.return = TRUE)$ix, n = generationSize)
  suggestions <- ga_res@population[best_idx, ]
  
  start_time <- Sys.time()
  ga_res <- ga(
    type = "real-valued",
    fitness = function(x) -obj_func(x, pars$names, 10),
    names = pars$names,
    lower = pars$lower,
    upper = pars$upper,
    popSize = generationSize,
    maxiter = generations,
    run = 5,
    parallel = cores,
    postFitness = function(obj, ...) postfit(obj, paste0(ga_iter, "_narrow")),
    monitor = TRUE,
    suggestions = suggestions,
    keepBest = TRUE
  )
  end_time <- Sys.time()
  print(end_time - start_time)
  
  save(
    ga_res,
    file = paste0(ga_final, "_narrow")
  )
  print(paste("Broad result saved on ", paste0(ga_final, "_narrow")))
}

# pQ in [0, 0.35] ####
pars <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ"),
  "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1, 0),
  "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5, 0.35)
)
ga_iter <- "data/several_runs/ga_range1_"
ga_final <- "data/several_runs/ga_range1_final.Rda"

run_ga(
  pars,
  ga_iter,
  ga_final
)

# pQ in [0.25, 0.60] ####
pars <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ"),
  "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1, 0.25),
  "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5, 0.6)
)
ga_iter <- "data/several_runs/ga_range2_"
ga_final <- "data/several_runs/ga_range2_final.Rda"

run_ga(
  pars,
  ga_iter,
  ga_final
)

# pQ in [0.5, 0.85] ####
pars <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ"),
  "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1, 0.5),
  "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5, 0.85)
)
ga_iter <- "data/several_runs/ga_range3_"
ga_final <- "data/several_runs/ga_range3_final.Rda"

run_ga(
  pars,
  ga_iter,
  ga_final
)

# pQ in [0.75, 1] ####
pars <- list(
  "names" = c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ"),
  "lower" = c(0.4,  0.025,  0.0001,  0.01,  3,  1,  1, 0.75),
  "upper" = c(0.9, 0.15, 0.03, 0.075, 12, 5, 5, 1)
)
ga_iter <- "data/several_runs/ga_range4_"
ga_final <- "data/several_runs/ga_range4_final.Rda"

run_ga(
  pars,
  ga_iter,
  ga_final
)

print("Done")

# Draft ####

load(paste0("data/several_runs/ga_range1_narrow_final.Rda"))
print(ga_res@solution)
load(paste0("data/one_cell_v2/ga_pQ_final.Rda"))
print(ga_res@solution)

get_best_sim_data <- function(sim_prefix, nr_sims = 1, base_path = base_path, func = run_one_cell) {
  load(paste0(base_path, "/", sim_prefix, "_final.Rda"))
  load("data/one_cell/weights.Rda")
  
  best_data <- func(
    ga_res@solution[1, ],
    return_result = TRUE,
    nr_sims = nr_sims,
    sim_parallel = TRUE
  )
  
  diversity_data <- bind_rows(
    Z13264 = do.call(
      rbind, 
      lapply(1:nr_sims, function(sim_idx) 
        get_diversity_df(best_data$Z13264[[sim_idx]]) %>%
          mutate(Sim_Idx = sim_idx)
      )
    ),
    Z14004 = do.call(
      rbind, 
      lapply(1:nr_sims, function(sim_idx) 
        get_diversity_df(best_data$Z14004[[sim_idx]]) %>%
          mutate(Sim_Idx = sim_idx)
      )
    ),
    .id = "Subject"
  )
  
  contribution_data <- bind_rows(
    Z13264 = do.call(
      rbind, 
      lapply(1:nr_sims, function(sim_idx) 
        get_contribution_dist_df(
          best_data$Z13264[[sim_idx]],
          break_list = break_list$Z13264
        ) %>%
          mutate(Sim_Idx = sim_idx)
      )
    ),
    Z14004 = do.call(
      rbind, 
      lapply(1:nr_sims, function(sim_idx) 
        get_contribution_dist_df(
          best_data$Z14004[[sim_idx]],
          break_list = break_list$Z14004
        ) %>%
          mutate(Sim_Idx = sim_idx)
      )
    ),
    .id = "Subject"
  )
  
  so_data <- bind_rows(
    Z13264 = do.call(
      rbind, 
      lapply(1:nr_sims, function(sim_idx) 
        get_so_df(best_data$Z13264[[sim_idx]]) %>%
          mutate(Sim_Idx = sim_idx)
      )
    ),
    Z14004 = do.call(
      rbind, 
      lapply(1:nr_sims, function(sim_idx) 
        get_so_df(best_data$Z14004[[sim_idx]]) %>%
          mutate(Sim_Idx = sim_idx)
      )
    ),
    .id = "Subject"
  )
  
  residuals_data <- do.call(
    rbind, 
    lapply(1:nr_sims, function(sim_idx) {
      df <- get_full_diffs_df(list(
        Z13264 = get_differences(list(best_data$Z13264[[sim_idx]]), 1, "Z13264"),
        Z14004 = get_differences(list(best_data$Z14004[[sim_idx]]), 1, "Z14004")
      )) %>%
        mutate(Sim_Idx = sim_idx)  %>%
        merge(
          point_weights,
          all.x = TRUE,
          all.y = FALSE
        ) %>%
        merge(
          graph_weights,
          by = c("Subject", "Graph"),
          all.x = TRUE,
          all.y = FALSE,
          suffixes = c("_Point", "_Graph")
        ) %>%
        mutate(
          Res = 
            Res*
            coalesce(Weight_Point, 1)*
            sqrt(coalesce(abs(Weight_Graph), 1))
        ) %>%
        mutate(Aspect = paste(Subject, Graph, sep = "_")) %>%
        group_by(Aspect) %>%
        mutate(RSS = sum(Res**2)) %>%
        ungroup()
      
      df$Total <- sum((df %>% group_by(Aspect) %>% summarise(RSS = RSS[1]))$RSS)
      
      return(df)
    })
  )
  
  return(list(
    Diversity = diversity_data,
    Contribution = contribution_data,
    SO = so_data,
    Residuals = residuals_data
  ))
}

nr_sims <- 50
range_1 <- get_best_sim_data("ga_range1_narrow", base_path = "data/several_runs", nr_sims = nr_sims)
range_2 <- get_best_sim_data("ga_range2_narrow", base_path = "data/several_runs", nr_sims = nr_sims)
range_3 <- get_best_sim_data("ga_range3_narrow", base_path = "data/several_runs", nr_sims = nr_sims)
range_4 <- get_best_sim_data("ga_range4_narrow", base_path = "data/several_runs", nr_sims = nr_sims)
best_two_full <- get_best_sim_data("ga_two", nr_sims = nr_sims, base_path = "data/one_cell")
best_two_reduced <- get_best_sim_data("ga_two", nr_sims = nr_sims, base_path = "data/one_cell_v2")
best_pQ_full <- get_best_sim_data("ga_pQ", nr_sims = nr_sims, base_path = "data/one_cell")
best_pQ_reduced <- get_best_sim_data("ga_pQ", nr_sims = nr_sims, base_path = "data/one_cell_v2")

save(
  range_1,
  range_2,
  range_3,
  range_4,
  best_two_full,
  best_two_reduced,
  best_pQ_full,
  best_pQ_reduced,
  
  file = "best_sims.Rda"
)

best_residuals <- bind_rows(
  "With pQ (06.01.25)" = best_pQ_full$Residuals,
  "With pQ (12.01.25)" = best_pQ_reduced$Residuals,
  "No pQ (06.01.25)" = best_two_full$Residuals,
  "No pQ (12.01.25)" = best_two_reduced$Residuals,
  # "No pA" = best_no_pA$Residuals,
  "Range 1" = range_1$Residuals,
  "Range 2" = range_2$Residuals,
  "Range 3" = range_3$Residuals,
  "Range 4" = range_4$Residuals,
  .id = "Source"
)

ggplot(
  best_residuals %>%
    group_by(Source, Total) %>%
    summarise(Source = Source[1], Total = Total[1]),
  aes(
    Source,
    Total
  )
) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point() +
  labs(
    title = "Total Residual Variance",
    y = "Total Residual"
  )

non_par_cols <- c(
  "Mean_Res",
  "Sd_Res",
  "First_Generation",
  "Last_Generation",
  "Total_Generations"
)
read_par_space <- function(
    sim_prefix, 
    base_path,
    min_gen = 1,
    max_gen = Inf
) {
  fl <- list.files(
    path = base_path,
    pattern = paste0(sim_prefix, "_[0-9]+\\.Rda")
  )
  
  do.call(rbind, lapply(fl, function(fn) {
    generation <- strsplit(fn, paste0(sim_prefix, "_"))[[1]][2]
    generation <- as.integer(strsplit(generation, ".Rda")[[1]][1])
    
    if(generation < min_gen | generation > max_gen) {
      return(NULL)
    }
    
    load(paste0(base_path, "/", fn))
    
    df <- data.frame(curr_ga@population)
    colnames(df) <- curr_ga@names
    df$Res <- -curr_ga@fitness
    df$Generation <- generation
    
    return(df)
  })) %>%
    group_by(across(c(-Res, - Generation))) %>% 
    summarise(
      Mean_Res = mean(Res),
      Sd_Res = sd(Res),
      First_Generation = min(Generation),
      Last_Generation = max(Generation),
      Total_Generations = length(Generation),
      .groups = "drop"
    )
}

par_space_1 <- read_par_space("ga_range1", base_path = "data/several_runs")
par_space_2 <- read_par_space("ga_range2", base_path = "data/several_runs")
par_space_3 <- read_par_space("ga_range3", base_path = "data/several_runs")
par_space_4 <- read_par_space("ga_range4", base_path = "data/several_runs")

ggplot(
  bind_rows(
    "Range 1" = par_space_1  %>%
      filter(Last_Generation == max(Last_Generation)),
    "Range 2" = par_space_2  %>%
      filter(Last_Generation == max(Last_Generation)),
    "Range 3" = par_space_3  %>%
      filter(Last_Generation == max(Last_Generation)),
    "Range 4" = par_space_4  %>%
      filter(Last_Generation == max(Last_Generation)),
    .id = "Model"
  ) %>%
    pivot_longer(!c(non_par_cols, "Model")),
  aes(
    x = value,
    color = paste0(Model, "(", Last_Generation, ")"),
    group = interaction(Model)
  )
) + 
  geom_density() +
  facet_wrap(vars(name), scales = "free") +
  labs(
    title = "Parameter Distributions for Last Generation",
    color = "Model (Final Generation)"
  )

GGally::ggpairs(
  bind_rows(
    "First_Gen" = bind_rows(
      par_space_1, 
      par_space_2, 
      par_space_3, 
      par_space_4
      ) %>% filter(
        First_Generation == 1
      ),
    "Range 1" = par_space_1 %>%
      slice_min(Mean_Res, n = 10),
    "Range 2" = par_space_2 %>%
      slice_min(Mean_Res, n = 10),
    "Range 3" = par_space_3 %>%
      slice_min(Mean_Res, n = 10),
    "Range 4" = par_space_4 %>%
      slice_min(Mean_Res, n = 10),
    .id = "Model"
  ),
  aes(
    color = Model,
    alpha = 0.1
  ),
  columns = 2:9,
  upper = "blank",
  diag = "blankDiag"
) + scale_color_manual(values = c(
  "First_Gen" = "lightgray",
  "Range 1" = "blue",
  "Range 2" = "red",
  "Range 3" = "green",
  "Range 4" = "yellow"
))
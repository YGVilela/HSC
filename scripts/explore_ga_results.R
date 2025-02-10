source("src/optimization.R")
library(ggplot2)
library(tidyverse)

base_path <- "data/one_cell_v2"

# Solution plots ####
get_best_sim_data <- function(sim_prefix, nr_sims = 1) {
  load(paste0(base_path, "/", sim_prefix, "_final.Rda"))
  load("data/one_cell/weights.Rda")
  
  best_data <- run_one_cell(
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

best_two <- get_best_sim_data("ga_two", nr_sims = 10)
best_pQ <- get_best_sim_data("ga_pQ", nr_sims = 10)

# Diversity plots ####
best_div <- bind_rows(
  "No pQ" = best_two$Diversity,
  "With pQ" = best_pQ$Diversity,
  .id = "Source"
)

ggplot(
  best_div %>% filter(Time > 365),
  aes(Time, Nr_Clones, color = Source, group = Source)
) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "min", geom = "line", linetype = "dotted") +
  stat_summary(fun = "max", geom = "line", linetype = "dotted") +
  geom_point(alpha = 0.2) +
  geom_point(
    data = clones_over_time %>% filter(Time > 365),
    mapping = aes(color = NULL, group = NULL),
    shape = "cross",
    size = 4,
    stroke = 1
  ) +
  facet_grid(cols = vars(Subject)) +
  scale_y_log10() +
  labs(
    title = "Diversity Fit",
    y = "Nr. Clones"
  )

# SO plots ####
best_so <- bind_rows(
  "No pQ" = best_two$SO,
  "With pQ" = best_pQ$SO,
  .id = "Source"
)

ggplot(
  best_so %>% filter(Total_Cells > 1000 & Time > 365),
  aes(Time, SO_Clone_Contribution, color = Source, group = Source)
) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "min", geom = "line", linetype = "dotted") +
  stat_summary(fun = "max", geom = "line", linetype = "dotted") +
  geom_point(alpha = 0.2) +
  geom_point(
    data = so_occurence %>% filter(Total_Cells > 1000 & Time > 365),
    mapping = aes(color = NULL, group = NULL),
    shape = "cross",
    size = 4,
    stroke = 1
  ) +
  facet_grid(cols = vars(Subject)) +
  scale_y_log10() +
  labs(
    title = "SO Clones Fit",
    y = "Fraction of SO Clones"
  )

# Contribution plots ####
best_cont <- bind_rows(
  "No pQ" = best_two$Contribution,
  "With pQ" = best_pQ$Contribution,
  .id = "Source"
)

ggplot(
  best_cont %>% filter(Time %in% c(362,397,462,495) & Clone_Contribution <= 1e-3),
  aes(Clone_Contribution, Frequency, color = Source, group = interaction(Source, Time))
) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "min", geom = "line", linetype = "dotted") +
  stat_summary(fun = "max", geom = "line", linetype = "dotted") +
  geom_point(alpha = 0.2) +
  geom_point(
    data = clone_contribution_dist %>% filter(Time %in% c(362,397,462,495) & Clone_Contribution <= 1e-3),
    mapping = aes(color = NULL, group = Time),
    shape = "cross",
    size = 4,
    stroke = 1
  ) +
  facet_wrap(vars(paste(Subject, Time, sep = "\n"))) +
  scale_y_log10() + scale_x_log10() +
  labs(
    title = "Clone Contribution Fit",
    x = "Clone Contribution"
  )

# Residuals ####
best_residuals <- bind_rows(
  "No pQ" = best_two$Residuals,
  "With pQ" = best_pQ$Residuals,
  .id = "Source"
)

# Final value ####
ggplot(
  best_residuals %>%
    group_by(Source, Total) %>%
    summarise(Source = Source[1], Total = Total[1]),
  aes(
    Source,
    Total
  )
) + 
  geom_point() +
  stat_summary(fun.data = function(dt) mean_sdl(dt, mult = 1)) +
  labs(
    title = "Total Residual Variance",
    y = "Total Residual"
  )

# Residual Distribution ####
ggplot(
  bind_rows(
    best_residuals,
    best_residuals %>% mutate(Aspect = "All")
  ),
  aes(
    Res,
    color = Aspect,
    group = interaction(Source, Aspect)
  )
) +
  geom_histogram(
    aes(y = after_stat(count)),
    binwidth = 0.03
  ) +
  geom_density(
    aes(y = after_stat(count)*0.03)
  ) +
  facet_grid(
    rows = vars(Aspect), 
    cols = vars(Source), 
    scales = "free_y"
  ) +
  labs(
    title = "Individual Residual Distribution",
    x = "Residual"
  )

# Single Point Residuals ####
ggplot(
  data = merge(
    best_residuals %>% filter(Source == "No pQ"),
    best_residuals %>% filter(Source == "With pQ"),
    by = c("Aspect", "X")
  ) %>% 
    group_by(Aspect, Res.x, Res.y) %>%
    summarise(
      Aspect = Aspect[1],
      Res.x = Res.x[1],
      Res.y = Res.y[1],
      .groups = "drop"
    ),
  aes(
    Res.x,
    Res.y,
    color = Aspect
  )
) +
  geom_abline(slope = 1, color = "red") +
  geom_point(alpha = 0.1) +
  geom_point(
    data =  merge(
      best_residuals %>% filter(Source == "No pQ") %>%
        group_by(Aspect, X, Source) %>%
        summarise(Res = mean(Res), .groups = "drop"),
      best_residuals %>% filter(Source == "With pQ") %>%
        group_by(Aspect, X, Source) %>%
        summarise(Res = mean(Res), .groups = "drop"),
      by = c("Aspect", "X")
    ),
    mapping = aes(
      Res.x, Res.y,
    ),
    shape = "cross",
    stroke = 2,
    color = "black"
  ) +
  facet_wrap(vars(Aspect)) +
  labs(
    title = "Single Point Residuals",
    x = "No pQ",
    y = "With pQ"
  )

# Aspect Residuals ####
ggplot(
  data = merge(
    best_residuals %>% filter(Source == "No pQ"),
    best_residuals %>% filter(Source == "With pQ"),
    by = c("Aspect")
  ) %>% 
    group_by(Aspect, RSS.x, RSS.y) %>%
    summarise(
      Aspect = Aspect[1],
      RSS.x = RSS.x[1],
      RSS.y = RSS.y[1],
      .groups = "drop"
    ),
  aes(
    RSS.x,
    RSS.y,
    color = Aspect
  )
) +
  geom_abline(slope = 1, color = "red")  +
  geom_point(alpha = 0.4) +
  geom_point(
    data = merge(
      best_residuals %>% filter(Source == "No pQ") %>%
        group_by(Aspect) %>%
        summarise(
          RSS = mean(RSS),
          .groups = "drop"
        ),
      best_residuals %>% filter(Source == "With pQ") %>%
        group_by(Aspect) %>%
        summarise(
          RSS = mean(RSS),
          .groups = "drop"
        ),
      by = c("Aspect")
    ),
    shape = "cross",
    stroke = 2,
    color = "black"
  ) +
  labs(
    title = "Aspect Residuals",
    x = "No pQ",
    y = "With pQ",
    color = "Aspect"
  )

# Param space plots ####
non_par_cols <- c(
  "Mean_Res",
  "Sd_Res",
  "First_Generation",
  "Last_Generation",
  "Total_Generations"
)
read_par_space <- function(sim_prefix, par_names) {
  fl <- list.files(
    path = base_path,
    pattern = paste0(sim_prefix, "_[0-9]+\\.Rda")
  )
  
  do.call(rbind, lapply(fl, function(fn) {
    generation <- strsplit(fn, paste0(sim_prefix, "_"))[[1]][2]
    generation <- as.integer(strsplit(generation, ".Rda")[[1]][1])
    
    load(paste0(base_path, "/", fn))
    
    df <- data.frame(pop)
    colnames(df) <- par_names
    df$Res <- -fit
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

par_space_two <- read_par_space(
  "ga_two",
  c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q")
)

par_space_pQ <- read_par_space(
  "ga_pQ",
  c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ")
)

# Param distribution for last generation ####
ggplot(
  bind_rows(
    "With pQ" = par_space_pQ,
    "No pQ" = par_space_two,
    .id = "Model"
  ) %>% 
    filter(Last_Generation == max(Last_Generation)) %>%
    pivot_longer(!c(non_par_cols, "Model")),
  aes(
    x = value,
    color = as.factor(Model),
    group = interaction(Model)
  )
) + 
  geom_density() +
  facet_wrap(vars(name), scales = "free") +
  labs(
    title = "Parameter Distributions for Last Generation",
    color = "Model"
  )

# Param distribution over time ####
ggplot(
  par_space_pQ %>% 
    filter(Mean_Res < 10) %>%
    filter(
      Last_Generation %in% c(
        min(Last_Generation),
        ceiling((max(Last_Generation) - min(Last_Generation))/2),
        max(Last_Generation)
      )
    ) %>%
    pivot_longer(!non_par_cols),
  aes(
    x = value,
    y = after_stat(density),
    color = as.factor(Last_Generation),
    group = interaction(Last_Generation)
  )
) + 
  geom_density() +
  facet_wrap(vars(name), scales = "free") +
  labs(
    title = "Parameter Distributions (With pQ)",
    color = "Generation"
  )

ggplot(
  par_space_two %>% 
    filter(Mean_Res < 10) %>%
    filter(
      Last_Generation %in% c(
        min(Last_Generation),
        ceiling((max(Last_Generation) - min(Last_Generation))/2),
        max(Last_Generation)
      )
    ) %>%
    pivot_longer(!Last_Generation & !First_Generation & !Total_Generations),
  aes(
    x = value,
    y = after_stat(density),
    color = as.factor(Last_Generation),
    group = interaction(Last_Generation)
  )
) + 
  geom_density() +
  facet_wrap(vars(name), scales = "free") +
  labs(
    title = "Parameter Distributions (No pQ)",
    color = "Generation"
  )

# Param pairs (With pQ) ####
{
  pairs_pQ <- do.call(rbind, expand.grid(
    colnames((par_space_pQ %>% select(!non_par_cols))),
    colnames((par_space_pQ %>% select(!non_par_cols)))
  ) %>% filter(Var1 != Var2) %>%
    apply(1, function(curr_row) {
      curr_data <- (par_space_pQ)[, c("Mean_Res", "Last_Generation", curr_row)]
      
      colnames(curr_data) <- c("Mean_Res", "Last_Generation", "Value1", "Value2")
      curr_data$Var1 <- curr_row["Var1"]
      curr_data$Var2 <- curr_row["Var2"]
      
      return(curr_data)
    }))
  
  relevant_percentiles <- c(0, 1, 5, 10, 50, 100)
  percentile_limits <- quantile(
    (par_space_pQ %>% filter(Last_Generation == max(Last_Generation)))$Mean_Res,
    relevant_percentiles/100
  )
  names(percentile_limits) <- c(
    paste0(relevant_percentiles[2:(length(relevant_percentiles) - 1)], "%"),
    paste0(">", relevant_percentiles[(length(relevant_percentiles) - 1)], "%"),
    "Shouldn't appear"
  )
  percentile_texts <- sapply(
    1:4, 
    function(idx) paste0(
      names(percentile_limits)[idx], ": [",
      round(percentile_limits[idx], digits = 2), ", ",
      round(percentile_limits[idx+1], digits = 2), ")"
    )
  )
}

# Residual distribution
{
  ggplot() +
    geom_density(
      data = par_space_pQ %>%
        filter(Last_Generation == max(Last_Generation)),
      mapping = aes(
        x = Mean_Res
      )
    ) +
    geom_vline(
      aes(
        xintercept = percentile_limits[2:5], 
        color = factor(percentile_texts, levels = percentile_texts)
      ),
      linetype = "dashed"
    ) +
    scale_x_log10() +
    labs(
      title = "Residual Distribution for Last Generation (With pQ)",
      x = "Residual",
      color = "Bounds:",
      y = "Count"
    )
}

# Full space
{
  ggplot(
    mapping = aes(Value1, Value2)
  ) +
    geom_point(
      data = pairs_pQ %>% filter(Last_Generation == 1),
      color = "lightgray"
    ) +
    geom_point(
      data = pairs_pQ %>% 
        filter(Last_Generation == max(Last_Generation)) %>%
        mutate(
          Pct_Idx = findInterval(
            Mean_Res,
            percentile_limits,
            rightmost.closed = TRUE
          ),
          Pct = names(percentile_limits)[Pct_Idx]
        ) %>%
        arrange(desc(Pct_Idx)),
      # color = "black"
      mapping = aes(
        color = factor(Pct, levels = names(percentile_limits))
      )
    ) +
    geom_point(
      data = pairs_pQ %>%
        filter(Mean_Res == min(Mean_Res)),
      color = "red",
      shape = "cross",
      stroke = 2
    ) +
    facet_grid(
      rows = vars(Var2),
      cols = vars(Var1),
      scales = "free"
    ) +
    scale_color_viridis_d() +
    labs(
      title = "Parameter Pairs on Last Generation (With pQ)",
      color = "Rank"
    )
}

# Param pairs (No pQ) ####
{
  pairs_two <- do.call(rbind, expand.grid(
    colnames((par_space_two %>% select(!non_par_cols))),
    colnames((par_space_two %>% select(!non_par_cols)))
  ) %>% filter(Var1 != Var2) %>%
    apply(1, function(curr_row) {
      curr_data <- (par_space_two)[, c("Mean_Res", "Last_Generation", curr_row)]
      
      colnames(curr_data) <- c("Mean_Res", "Last_Generation", "Value1", "Value2")
      curr_data$Var1 <- curr_row["Var1"]
      curr_data$Var2 <- curr_row["Var2"]
      
      return(curr_data)
    }))
  
  relevant_percentiles <- c(0, 1, 5, 10, 50, 100)
  percentile_limits <- quantile(
    (par_space_two %>% filter(Last_Generation == max(Last_Generation)))$Mean_Res,
    relevant_percentiles/100
  )
  names(percentile_limits) <- c(
    paste0(relevant_percentiles[2:(length(relevant_percentiles) - 1)], "%"),
    paste0(">", relevant_percentiles[(length(relevant_percentiles) - 1)], "%"),
    "Shouldn't appear"
  )
  percentile_texts <- sapply(
    1:4, 
    function(idx) paste0(
      names(percentile_limits)[idx], ": [",
      round(percentile_limits[idx], digits = 2), ", ",
      round(percentile_limits[idx+1], digits = 2), ")"
    )
  )
}

# Residual distribution
{
  ggplot() +
    geom_density(
      data = par_space_two %>%
        filter(Last_Generation == max(Last_Generation)),
      mapping = aes(
        x = Mean_Res
      )
    ) +
    geom_vline(
      aes(
        xintercept = percentile_limits[2:5], 
        color = factor(percentile_texts, levels = percentile_texts)
      ),
      linetype = "dashed"
    ) +
    scale_x_log10() +
    labs(
      title = "Residual Distribution for Last Generation (No pQ)",
      x = "Residual",
      color = "Bounds:",
      y = "Count"
    )
}

# Full space
{
  ggplot(
    mapping = aes(Value1, Value2)
  ) +
    geom_point(
      data = pairs_two %>% filter(Last_Generation == 1),
      color = "lightgray"
    ) +
    geom_point(
      data = pairs_two %>% 
        filter(Last_Generation == max(Last_Generation)) %>%
        mutate(
          Pct_Idx = findInterval(
            Mean_Res,
            percentile_limits,
            rightmost.closed = TRUE
          ),
          Pct = names(percentile_limits)[Pct_Idx]
        ) %>%
        arrange(desc(Pct_Idx)),
      # color = "black"
      mapping = aes(
        color = factor(Pct, levels = names(percentile_limits))
      )
    ) +
    geom_point(
      data = pairs_two %>%
        filter(Mean_Res == min(Mean_Res)),
      color = "red",
      shape = "cross",
      stroke = 2
    ) +
    facet_grid(
      rows = vars(Var2),
      cols = vars(Var1),
      scales = "free"
    ) +
    scale_color_viridis_d() +
    labs(
      title = "Parameter Pairs on Last Generation (No pQ)",
      color = "Rank"
    )
}
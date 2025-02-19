source("src/optimization.R")
library(ggplot2)
library(tidyverse)

base_path <- "data/one_cell_v2"

# Solution plots ####
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

st <- Sys.time()
best_two_full <- get_best_sim_data("ga_two", nr_sims = 10, base_path = "data/one_cell")
best_two_reduced <- get_best_sim_data("ga_two", nr_sims = 10, base_path = "data/one_cell_v2")
best_pQ_full <- get_best_sim_data("ga_pQ", nr_sims = 10, base_path = "data/one_cell")
best_pQ_reduced <- get_best_sim_data("ga_pQ", nr_sims = 10, base_path = "data/one_cell_v2")
best_no_pA <- get_best_sim_data("ga", nr_sims = 10, base_path = "data/no_pA", func = run_no_pA)
et <- Sys.time()
print(et - st)

# save(best_two, best_pQ, file = "something.Rda")

# Diversity plots ####
best_div <- bind_rows(
  "With pQ (06.01.25)" = best_pQ_full$Diversity,
  "With pQ (12.01.25)" = best_pQ_reduced$Diversity,
  "No pQ (06.01.25)" = best_two_full$Diversity,
  "No pQ (12.01.25)" = best_two_reduced$Diversity,
  "No pA" = best_no_pA$Diversity,
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
  "With pQ (06.01.25)" = best_pQ_full$SO,
  "With pQ (12.01.25)" = best_pQ_reduced$SO,
  "No pQ (06.01.25)" = best_two_full$SO,
  "No pQ (12.01.25)" = best_two_reduced$SO,
  "No pA" = best_no_pA$SO,
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
  "With pQ (06.01.25)" = best_pQ_full$Contribution,
  "With pQ (12.01.25)" = best_pQ_reduced$Contribution,
  "No pQ (06.01.25)" = best_two_full$Contribution,
  "No pQ (12.01.25)" = best_two_reduced$Contribution,
  "No pA" = best_no_pA$Contribution,
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
  "With pQ (06.01.25)" = best_pQ_full$Residuals,
  "With pQ (12.01.25)" = best_pQ_reduced$Residuals,
  "No pQ (06.01.25)" = best_two_full$Residuals,
  "No pQ (12.01.25)" = best_two_reduced$Residuals,
  "No pA" = best_no_pA$Residuals,
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
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_point() +
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
best_residuals %>%
  group_by(Aspect, Source, X) %>%
  summarise(Res = mean(Res), .groups = "drop") %>%
  pivot_wider(names_from = Source, values_from = Res) %>%
  GGally::ggpairs(
    mapping = aes(colour = Aspect),
    diag = list(continuous = function(data, mapping, ...) {
      mapping$colour <- NULL
      GGally::ggally_densityDiag(data, mapping)
    }),
    lower = list(continuous = function(data, mapping, ...) {
      GGally::ggally_points(data, mapping) +
        geom_abline(slope = 1, color = "red")
    }),
    upper = list(continuous = "blank"), 
    columns = 3:10,
    legend = c(2, 1)
  )

best_residuals %>%
  group_by(Aspect, Source) %>%
  summarise(RSS = mean(RSS), .groups = "drop") %>%
  pivot_wider(names_from = Source, values_from = RSS) %>%
  GGally::ggpairs(
    mapping = aes(colour = Aspect),
    diag = list(continuous = "blankDiag"),
    lower = list(continuous = function(data, mapping, ...) {
      GGally::ggally_points(data, mapping) +
        geom_abline(slope = 1, color = "red")
    }),
    upper = list(continuous = "blank"), 
    columns = 2:9,
    legend = c(2, 1)
  )

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

read_par_space_old <- function(
    sim_prefix,
    par_names,
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

par_space_no_pA <- read_par_space(
  "ga",
  base_path = "data/no_pQ"
)

par_space_two <- read_par_space(
  "ga_two",
  base_path = "data/one_cell_v2"
)

par_space_pQ <- read_par_space(
  "ga_pQ",
  base_path = "data/one_cell_v2"
)

par_space_pQ_old <- read_par_space_old(
  "ga_pQ",
  c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q", "pQ"),
  base_path = "data/one_cell",
  max_gen = 10
)

par_space_two_old <- read_par_space_old(
  "ga_two",
  c("pA", "dA", "tQA", "tAQ", "clone_mult", "space_mult_A", "space_mult_Q"),
  base_path = "data/one_cell",
  max_gen = 50
)

par_space_no_pA <- read_par_space("ga", base_path = "data/no_pQ")

# Param distribution for last generation ####
ggplot(
  bind_rows(
    "With pQ" = par_space_pQ %>%
      # filter(Last_Generation == max(Last_Generation)),
      slice_min(Mean_Res, n = 100),
    "No pQ" = par_space_two %>%
      # filter(Last_Generation == max(Last_Generation)),
      slice_min(Mean_Res, n = 100),

    "With pQ (old)" = par_space_pQ_old %>%
      # filter(Last_Generation == max(Last_Generation)),
      slice_min(Mean_Res, n = 100),
    "No pQ (old)" = par_space_two_old %>%
      # filter(Last_Generation == max(Last_Generation)),
      slice_min(Mean_Res, n = 100),
    .id = "Model"
  ) %>%
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
# With pQ
{
  ggplot(
    par_space_pQ_old %>% 
      filter(Mean_Res < 10) %>%
      filter(
        Last_Generation == 10
      ) %>%
      pivot_longer(!non_par_cols),
    aes(
      x = value,
      y = after_stat(density)
      # color = as.factor(Last_Generation),
      # group = interaction(Last_Generation)
    )
  ) + 
    geom_density() +
    geom_density(data = par_space_pQ_old %>% 
                   filter(Mean_Res < 10) %>%
                   filter(
                     First_Generation == 1
                   ) %>%
                   pivot_longer(!non_par_cols), color = "red") +
    facet_wrap(vars(name), scales = "free") +
    labs(
      title = "Parameter Distributions (With pQ)",
      color = "Generation"
    )
}

# No pQ
{
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
      title = "Parameter Distributions (No pQ)",
      color = "Generation"
    )
}

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
    (par_space_pQ %>%
       filter(Last_Generation >= max(Last_Generation) - 2))$Mean_Res,
       # filter(Last_Generation == max(Last_Generation)))$Mean_Res,
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
        # filter(Last_Generation == max(Last_Generation)),
        filter(Last_Generation >= max(Last_Generation) - 2),
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
        # filter(Last_Generation == max(Last_Generation)) %>%
        filter(Last_Generation >= max(Last_Generation) - 2) %>%
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
        # filter(Last_Generation >= max(Last_Generation) - 0) %>%
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
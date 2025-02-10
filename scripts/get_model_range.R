source("src/optimization.R")

sim_results_file <- "data/one_cell/par_space_diffs.Rda"

if(!file.exists(sim_results_file)) {
  nr_samples <- 1000
  cores <- min(100, parallel::detectCores())
  
  # Two compartments (simulations) ####
  sample_pars <- data.frame(
    pA = sample(seq(0.4, 0.9, length.out = nr_samples)),
    dA = sample(seq(0.025, 0.15, length.out = nr_samples)),
    tAQ = sample(seq(0.01, 0.075, length.out = nr_samples)),
    tQA = sample(seq(0.0001, 0.03, length.out = nr_samples)),
    clone_mult = sample(seq(3, 12, length.out = nr_samples)),
    space_mult_A = sample(seq(1, 5, length.out = nr_samples)),
    space_mult_Q = sample(seq(1, 5, length.out = nr_samples))
  )
  
  start_time <- Sys.time()
  diffs <- parallel::mclapply(
    1:nr_samples, 
    function(idx) 
      tryCatch(
        {
          print(paste("Starting", idx))
          start_time <- Sys.time()
          
          ret <- run_one_cell(as.list(sample_pars[idx, ]))
          
          end_time <- Sys.time()
          print(paste("Done with", idx))
          print(end_time - start_time)
          
          return(ret)
        },
        error = function(e) {
          print(paste("Error in", idx))
          print(conditionMessage(e))
          
          return(NULL)
        } 
      ),
    mc.cores = cores
  )
  
  end_time <- Sys.time()
  print("Done with sims")
  print(end_time - start_time)
  
  full_df <- do.call(rbind, lapply(1:nr_samples, function(idx) {
    tryCatch(
      {
        if(is.null(diffs[[idx]])) {
          print(paste(idx, "is NULL"))
          return(NULL)
        }
        
        df <- get_full_diffs_df(diffs[[idx]])
        
        df$Sim_Idx <- idx
        
        return(df)
      },
      error = function(e) {
        print(paste("Error in", idx))
        print(conditionMessage(e))
        return(NULL)
      }
    )
    
  }))
  
  save(
    sample_pars, full_df,
    file = sim_results_file
  )
} else {
  print("Simulation data already exists.")
}

stop("Ok")

# Plotting results ####

library(ggplot2)
library(gghalves)
library(dplyr)

get_all_plots <- function(file_name, plt_title, run_func = NULL, top_number = 100) {
  load(file_name)
  
  # Residual Distribution ####
  adjusted_full_df <- full_df %>%
    group_by(Subject, Graph, X) %>%
    mutate(
      SD_Res = sd(Res, na.rm = TRUE),
      Mean_Res = mean(Res, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      Normalized_Res = (Res)/(ifelse(SD_Res > 0, 2*SD_Res, 1))
    )
  
  adjusted_full_df[grepl("Contribution", adjusted_full_df$Graph), "X"] <- 
    log10(adjusted_full_df[grepl("Contribution", adjusted_full_df$Graph), "X"])
  
  residual_plot <- ggplot(
    adjusted_full_df,
    aes(
      x = as.factor(X),
      y = Res,
      group = interaction(X, Graph, Subject)
    )
  ) + 
    stat_summary(fun.data = "mean_sdl", geom = "errorbar") +
    geom_half_violin(side = "r", draw_quantiles = c(0.25, 0.5, 0.75), linetype = "dashed") +
    stat_summary(fun = "mean", geom = "point", color = "red") +
    facet_wrap(
      vars(interaction(Subject, Graph, sep = "\n")),
      scales = "free"
    ) + labs(
      title = paste("Residual Distribution -", plt_title),
      x = "X",
      y = "Residual"
    )
  
  normalized_residual_plot <- ggplot(
    adjusted_full_df,
    aes(
      x = as.factor(X),
      y = Normalized_Res,
      group = interaction(X, Graph, Subject)
    )
  ) + 
    stat_summary(fun.data = "mean_sdl", geom = "errorbar") +
    geom_half_violin(side = "r", draw_quantiles = c(0.25, 0.5, 0.75), linetype = "dashed") +
    stat_summary(fun = "mean", geom = "point", color = "red") +
    facet_wrap(
      vars(interaction(Subject, Graph, sep = "\n")),
      scales = "free"
    ) + labs(
      title = paste("Residual Distribution (Normalized) -", plt_title),
      x = "X",
      y = "Residual"
    )
  
  # RSS Distribution ####
  rss_df <- adjusted_full_df %>%
    group_by(Subject, Graph, Sim_Idx) %>%
    summarise(RSS = sum(Normalized_Res**2)) %>%
    group_by(Subject, Graph) %>%
    mutate(
      SD_RSS = sd(RSS, na.rm = TRUE),
      Mean_RSS = mean(RSS, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      Normalized_RSS = (RSS)/(2*SD_RSS)
    ) %>%
    ungroup()
  
  rss_plot <- ggplot(
    rss_df,
    aes(
      x = interaction(Subject, Graph, sep = "\n"),
      y = RSS
    )
  ) + 
    stat_summary(fun.data = "mean_sdl", geom = "errorbar") +
    geom_half_violin(side = "r", draw_quantiles = c(0.25, 0.5, 0.75), linetype = "dashed") +
    stat_summary(fun = "mean", geom = "point", color = "red")  +
    labs(
      title = paste("RSS Distribution-", plt_title),
      x = "Graph",
      y = "RSS"
    )
  
  normalized_rss_plot <- ggplot(
    rss_df,
    aes(
      x = interaction(Subject, Graph, sep = "\n"),
      y = Normalized_RSS
    )
  ) + 
    stat_summary(fun.data = "mean_sdl", geom = "errorbar") +
    geom_half_violin(side = "r", draw_quantiles = c(0.25, 0.5, 0.75), linetype = "dashed") +
    stat_summary(fun = "mean", geom = "point", color = "red")  +
    labs(
      title = paste("RSS Distribution (Normalized) -", plt_title),
      x = "Graph",
      y = "RSS"
    )
  
  # Final Sum Distribution ####
  final_value_df <- rss_df %>%
    group_by(Sim_Idx) %>%
    summarise(Final_Value = sum(Normalized_RSS)) %>%
    ungroup()
  
  final_sum_plot <- ggplot(
    data = final_value_df,
    aes(
      x = Final_Value
    )
  ) + geom_density() +
    labs(
      title = paste("Final Sum Distribution -", plt_title),
      x = "Final Sum",
      y = "Density"
    )
  
  # Par space points ####
  pars_with_residuals <- cbind(sample_pars[final_value_df$Sim_Idx,], Res = final_value_df$Final_Value) %>%
    arrange(Res)
  parameter_pairs <- do.call(rbind, expand.grid(
    colnames((pars_with_residuals %>% select(!Res))),
    colnames((pars_with_residuals %>% select(!Res)))
  ) %>% filter(Var1 != Var2) %>%
    apply(1, function(curr_row) {
      curr_data <- (pars_with_residuals)[, c("Res", curr_row)]
      
      colnames(curr_data) <- c("Res", "Value1", "Value2")
      curr_data$Var1 <- curr_row["Var1"]
      curr_data$Var2 <- curr_row["Var2"]
      
      return(curr_data)
    }))
  
  
  best_pars <- pars_with_residuals %>% 
    slice_min(Res, n = top_number) %>%
    mutate(Fit  = 1:length(Res))
  
  best_pars_pairs <- do.call(rbind, expand.grid(
    colnames((best_pars %>% select(!Res & !Fit))),
    colnames((best_pars %>% select(!Res & !Fit)))
  ) %>% filter(Var1 != Var2) %>%
    apply(1, function(curr_row) {
      curr_data <- (best_pars)[, c("Res", "Fit", curr_row)]
      
      colnames(curr_data) <- c("Res", "Fit", "Value1", "Value2")
      curr_data$Var1 <- curr_row["Var1"]
      curr_data$Var2 <- curr_row["Var2"]
      
      return(curr_data)
    }))
  
  par_space_plot <- ggplot(
    mapping = aes(
      x = Value1,
      y = Value2,
      color = Res
    )
  ) +
    geom_point(
      data = parameter_pairs %>% filter(
        Res - min(Res) > 0.5*(max(Res) - min(Res))
      ),
      color = "lightgray"
    ) +
    geom_point(
      data = parameter_pairs %>% filter(
        Res - min(Res) <= 0.5*(max(Res) - min(Res))
      ),
      color = "black"
    ) +
    geom_point(
      data = best_pars_pairs,
      mapping = aes(
        color = paste(strsplit(intToUtf8(64 + Fit), "")[[1]], round(Res, digits = 2), sep = ": ")
      ),
      shape = "cross",
      stroke = 2
    ) +
    facet_grid(
      cols = vars(Var1),
      rows = vars(Var2),
      scales = "free"
    ) + labs(
      color = "Best Plots",
      title = paste("Parameter Space -", plt_title),
    )
  
  # Best plots ####
  if(!is.null(run_func)) {
    best_data <- parallel::mclapply(
      1:nrow(best_pars), 
      function(idx) res <- run_func(best_pars[idx, ]),
      mc.cores = parallel::detectCores()
    )
    
    best_df <- do.call(rbind, lapply(1:length(best_data), function(idx) {
      df <- bind_rows(
        #### Z13264 ####
        Z13264 = bind_rows(
          Diversity = best_data[[idx]]$Z13264$Diversity %>% rename(
            X = Time,
            Y = Nr_Clones_Sim,
            Y_Data = Nr_Clones_Data
          ),
          SO = best_data[[idx]]$Z13264$SO %>% rename(
            X = Time,
            Y = SO_Clone_Contribution_Sim,
            Y_Data = SO_Clone_Contribution_Data
          ),
          Contribution_397 = best_data[[idx]]$Z13264$Contribution %>%
            filter(Time == 397) %>%
            rename(
              X = Clone_Contribution,
              Y = Frequency_Sim,
              Y_Data = Frequency_Data
            ),
          Contribution_462 = best_data[[idx]]$Z13264$Contribution %>%
            filter(Time == 462) %>%
            rename(
              X = Clone_Contribution,
              Y = Frequency_Sim,
              Y_Data = Frequency_Data
            ),
          .id = "Graph"
        ),
        
        #### Z14004 ####
        Z14004 = bind_rows(
          Diversity = best_data[[idx]]$Z14004$Diversity %>% rename(
            X = Time,
            Y = Nr_Clones_Sim,
            Y_Data = Nr_Clones_Data
          ),
          SO = best_data[[idx]]$Z14004$SO %>% rename(
            X = Time,
            Y = SO_Clone_Contribution_Sim,
            Y_Data = SO_Clone_Contribution_Data
          ),
          Contribution_362 = best_data[[idx]]$Z14004$Contribution %>%
            filter(Time == 362) %>%
            rename(
              X = Clone_Contribution,
              Y = Frequency_Sim,
              Y_Data = Frequency_Data
            ),
          Contribution_495 = best_data[[idx]]$Z14004$Contribution %>%
            filter(Time == 495) %>%
            rename(
              X = Clone_Contribution,
              Y = Frequency_Sim,
              Y_Data = Frequency_Data
            ),
          
          # End ####
          .id = "Graph"
        ),
        
        .id = "Subject"
      ) %>% select(
        Subject, Graph,
        X, Y, Res, Y_Data
      )
      
      df$Sim_Idx <- idx
      
      return(df)
    }))
    
    best_graph_plots <- ggplot(
      best_df,
      aes(
        x = X,
        y = Y,
        color = strsplit(intToUtf8(64 + Sim_Idx), "")[[1]],
        linetype = "Simulation",
        group = interaction(Sim_Idx, Subject)
      )
    ) +
      geom_line() +
      geom_point() +
      geom_line(
        mapping = aes(
          y = Y_Data,
          linetype = "Data",
          color = NULL,
          group = Subject
        ),
        color = "black"
      ) +
      facet_wrap(
        vars(paste0(Subject, "\n", Graph)),
        scales = "free",
        nrow = 2
      ) +
      scale_linetype_manual(values = c(
        "Data" = "dashed",
        "Simulation" = "solid"
      )) +
      scale_x_log10() + scale_y_log10()  + labs(
        color = "Best Plots",
        linetype = "Source",
        title = paste("Best Plots -", plt_title),
      )
  } else{
    best_graph_plots = NULL
  }
  
  
  # Return Graphs ####
  return(list(
    Residual = residual_plot,
    Normalized_Residual = normalized_residual_plot,
    RSS = rss_plot,
    Normalized_RSS = normalized_rss_plot,
    Final_Sum = final_sum_plot,
    Par_Space = par_space_plot,
    Best_Graphs = best_graph_plots
  ))
}

get_normalized_residuals <- function(file_name) {
  load(file_name)
  
  # Residual Distribution ####
  adjusted_full_df <- full_df %>%
    group_by(Subject, Graph, X) %>%
    mutate(
      SD_Res = sd(Res, na.rm = TRUE),
      Mean_Res = mean(Res, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      Normalized_Res = (Res)/(ifelse(SD_Res > 0, 2*SD_Res, 1))
    )
  
  rss_df <- adjusted_full_df %>%
    group_by(Subject, Graph, Sim_Idx) %>%
    summarise(RSS = sum(Normalized_Res**2)) %>%
    group_by(Subject, Graph) %>%
    mutate(
      SD_RSS = sd(RSS, na.rm = TRUE),
      Mean_RSS = mean(RSS, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      Normalized_RSS = (RSS)/(2*SD_RSS)
    ) %>%
    ungroup()
  
  fully_normalized_residues <- 
    merge(adjusted_full_df, rss_df) %>%
    mutate(Final_Point_Res = Normalized_Res/sqrt(2*SD_RSS)) %>%
    group_by(Sim_Idx) %>%
    mutate(
      Total_Res = sum(Final_Point_Res**2)
    ) %>%
    # select(
    #   Sim_Idx, Total_Res,
    #   Subject, Graph, X, Final_Point_Res
    # ) %>%
    ungroup()
  
  return(fully_normalized_residues)
}

two_compartment_plots <- get_all_plots(
  sim_results_file,
  "Two Compartments",
  run_func = run_one_cell,
  top_number = 5
)

save(
  two_compartment_plots,
  file = "data/one_cell/residual_weights_plots.Rda"
)

get_point_weights <- function(file_name) {
  load(file_name)
  
  # Residual Distribution ####
  adjusted_full_df <- full_df %>%
    group_by(Subject, Graph, X) %>%
    mutate(
      SD_Res = sd(Res, na.rm = TRUE),
      Mean_Res = mean(Res, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      Normalized_Res = (Res)/(ifelse(SD_Res > 0, 2*SD_Res, 1)),
      Weight = 1/ifelse(SD_Res > 0, 2*SD_Res, 1)
    )
  
  # RSS Distribution ####
  rss_df <- adjusted_full_df %>%
    group_by(Subject, Graph, Sim_Idx) %>%
    summarise(RSS = sum(Normalized_Res**2)) %>%
    group_by(Subject, Graph) %>%
    mutate(
      SD_RSS = sd(RSS, na.rm = TRUE),
      Mean_RSS = mean(RSS, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    mutate(
      Normalized_RSS = (RSS)/(2*SD_RSS),
      Weight = 1/(2*SD_RSS)
    ) %>%
    ungroup()
  
  # Return weights ####
  
  return(list(
    Point_Weight = adjusted_full_df %>%
      group_by(Subject, Graph, X) %>%
      summarise(Weight = Weight[1]) %>%
      select(Subject, Graph, X, Weight) %>%
      ungroup(),
    Graph_Weight = rss_df %>%
      group_by(Subject, Graph) %>%
      summarise(Weight = Weight[1]) %>%
      select(Subject, Graph, Weight) %>%
      ungroup()
  ))
}

two_point_weights <- get_point_weights(sim_results_file)
point_weights <- two_point_weights$Point_Weight
graph_weights <- two_point_weights$Graph_Weight

save(
  point_weights, graph_weights,
  file = "data/one_cell/weights.Rda"
)
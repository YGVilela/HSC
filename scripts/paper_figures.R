# Colors ####
# Base: https://coolors.co/e9727c-024150-fcd0a1-b1b695-53917e

color_scheme <- list(
  Subject = setNames( 
    c(
      '#9F1924', '#036277'
    ),
    c(
      "Z13264", "Z14004"
    )
  ),
  Aspects = setNames( 
    c(
      '#FCD0A1', '#B1B695', '#93987C','#53917E'
    ),
    c(
      "Clone Diversity", "Clone Size Dist. #1", "Clone Size Dist. #2", "SO Clone Fraction"
    )
  )
)

# Src ####
library(ggplot2)
library(ggpattern)
library(dplyr)
library(tidyr)
library(RColorBrewer)
source("src/Simplified_HSC_Model.R")

load("data/derived_radtke.Rda")
load("data/extracted_radtke.Rda")


# Metric plots ####
get_div_plot <- function(
    data_list,
    facet_by_model = TRUE,
    show_data = TRUE,
    add_pattern = TRUE,
    crop_y = TRUE
) {
  model_data <- bind_rows(
    data_list,
    .id = "Model"
  )
  
  div_plot <- ggplot(
    model_data,
    aes(
      Time, Nr_Clones,
      color = Subject,
      group = interaction(Subject, Model)
    )
  ) +
    stat_summary(fun = "mean", geom = "line") +
    stat_summary(fun = "mean", geom = "point") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_color_manual(
      values = color_scheme$Subject
    ) +
    labs(
      title = "Clonal Diversity",
      y = "Number of Clones",
      x = "Time [days]"
    )
  
  if(show_data) {
    div_plot <- div_plot + 
      geom_point(
        data = clones_over_time,
        aes(group = Subject),
        shape = "cross",
        stroke = 1.1
      ) +
      geom_line(
        data = clones_over_time,
        aes(group = Subject),
        linetype = "dashed"
      )
  }
  
  if(add_pattern) {
    div_plot <- div_plot +
      geom_rect_pattern(
        data = data.frame(c(1)),
        aes(
          xmin = 0,
          ymin = 0,
          xmax = 360,
          ymax = 20000
        ),
        inherit.aes = FALSE,
        pattern = 'stripe',
        pattern_angle = 45,
        pattern_density = 0.5,
        pattern_fill = 'white',
        fill = 'lightgray',
        pattern_alpha = 0.2,
        alpha = 0.2
      )
  }
  
  if(facet_by_model) {
    div_plot <- div_plot + facet_grid(cols = vars(Model))
  }
  
  if(crop_y) {
    div_plot <- div_plot + ylim(0, 20000)
  }
  
  return(div_plot)
}

get_so_plot <- function(
    data_list,
    facet_by_model = TRUE,
    show_data = TRUE,
    add_pattern = TRUE,
    crop_y = TRUE
) {
  model_data <- bind_rows(
    data_list,
    .id = "Model"
  ) %>%
    filter(Total_Cells > 100)
  
  so_plot <- ggplot(
    model_data,
    aes(
      Time, SO_Clone_Contribution,
      color = Subject,
      group = interaction(Subject, Model)
    )
  ) +
    stat_summary(fun = "mean", geom = "line") +
    stat_summary(fun = "mean", geom = "point") +
    theme_classic() +
    theme(legend.position = "none") +
    scale_color_manual(
      values = color_scheme$Subject
    ) +
    labs(
      title = "Single-time Occurring Clones",
      y = "Fraction of SO Clones",
      x = "Time [days]"
    )
  
  if(show_data) {
    so_plot <- so_plot + 
      geom_point(
        data = so_occurence %>% filter(Total_Cells > 100),
        aes(group = Subject),
        shape = "cross",
        stroke = 1.1
      ) +
      geom_line(
        data = so_occurence %>% filter(Total_Cells > 100),
        aes(group = Subject),
        linetype = "dashed"
      )
  }
  
  if(add_pattern) {
    so_plot <- so_plot +
      geom_rect_pattern(
        data = data.frame(c(1)),
        aes(
          xmin = 0,
          ymin = -0.01,
          xmax = 360,
          ymax = 0.75,
        ),
        inherit.aes = FALSE,
        pattern = 'stripe',
        pattern_angle = 45,
        pattern_density = 0.5,
        pattern_fill = 'white',
        fill = 'lightgray',
        pattern_alpha = 0.3,
        alpha = 0.3
      )
  }
  
  if(facet_by_model) {
    so_plot <- so_plot + facet_grid(cols = vars(Model))
  }
  
  if(crop_y) {
    so_plot <- so_plot + ylim(-0.01, 0.75)
  }
  
  return(so_plot)
}

get_contribution_Z14004_plot <- function(
    data_list,
    facet_by_model = TRUE,
    show_data = TRUE,
    add_pattern = TRUE,
    time_colors = c(
      "362" = "#04839F",
      "495" = "#024150"
    )
) {
  model_data <- bind_rows(
    data_list,
    .id = "Model"
  ) %>%
    filter(
      Time %in% as.integer(names(time_colors)) &
      Subject == "Z14004"
    )
  
  cont_plot <- ggplot(
    model_data,
    aes(
      Clone_Contribution, Frequency,
      color = as.factor(Time),
      group = interaction(Time, Model)
    )
  ) +
    stat_summary(fun = "mean", geom = "line") +
    stat_summary(fun = "mean", geom = "point")  +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_log10() + scale_x_log10() +
    scale_color_manual(
      values = time_colors
    ) +
    labs(
      title = "Clone Size Distribution (Z14004)",
      y = "Frequency",
      x = "Relative Clone Size"
    )
  
  if(show_data) {
    cont_plot <- cont_plot + 
      geom_point(
        data = clone_contribution_dist  %>% filter(
          Time %in% as.integer(names(time_colors)) &
          Subject == "Z14004"
        ),
        mapping = aes(group = Time),
        shape = "cross",
        stroke = 1.1
      ) +
      geom_line(
        data = clone_contribution_dist  %>% filter(
          Time %in% as.integer(names(time_colors)) &
          Subject == "Z14004"
        ),
        mapping = aes(group = Time),
        linetype = "dashed"
      )
  }
  
  if(add_pattern) {
    cont_plot <- cont_plot +
      geom_rect_pattern(
        data = data.frame(c(1)),
        aes(
          xmin = 1e-3,
          ymin = 0,
          xmax = 0.03981072,
          ymax = 0.3111111,
        ),
        inherit.aes = FALSE,
        pattern = 'stripe',
        pattern_angle = 45,
        pattern_density = 0.5,
        pattern_fill = 'white',
        fill = 'lightgray',
        pattern_alpha = 0.3,
        alpha = 0.3
      )
  }
  
  if(facet_by_model) {
    cont_plot <- cont_plot + facet_grid(cols = vars(Model))
  }
  
  return(cont_plot)
}

get_contribution_Z13264_plot <- function(
    data_list,
    facet_by_model = TRUE,
    show_data = TRUE,
    add_pattern = TRUE,
    time_colors = c(
      "397" = "#C21E2C",
      "462" = "#7C131B"
    )
) {
  model_data <- bind_rows(
    data_list,
    .id = "Model"
  ) %>%
    filter(
      Time %in% as.integer(names(time_colors)) &
      Subject == "Z13264"
    )
  
  cont_plot <- ggplot(
    model_data,
    aes(
      Clone_Contribution, Frequency,
      color = as.factor(Time),
      group = interaction(Time, Model)
    )
  ) +
    stat_summary(fun = "mean", geom = "line") +
    stat_summary(fun = "mean", geom = "point")  +
    theme_classic() +
    theme(legend.position = "none") +
    scale_y_log10() + scale_x_log10() +
    scale_color_manual(
      values = time_colors
    ) +
    labs(
      title = "Clone Size Distribution (Z13264)",
      y = "Frequency",
      x = "Relative Clone Size"
    )
  
  if(show_data) {
    cont_plot <- cont_plot + 
      geom_point(
        data = clone_contribution_dist  %>% filter(
          Time %in% as.integer(names(time_colors)) &
          Subject == "Z13264"
        ),
        mapping = aes(group = Time),
        shape = "cross",
        stroke = 1.1
      ) +
      geom_line(
        data = clone_contribution_dist  %>% filter(
          Time %in% as.integer(names(time_colors)) &
          Subject == "Z13264"
        ),
        mapping = aes(group = Time),
        linetype = "dashed"
      )
  }
  
  if(add_pattern) {
    cont_plot <- cont_plot +
      geom_rect_pattern(
        data = data.frame(c(1)),
        aes(
          xmin = 1e-3,
          ymin = 0,
          xmax = 0.01,
          ymax = 0.2694335,
        ),
        inherit.aes = FALSE,
        pattern = 'stripe',
        pattern_angle = 45,
        pattern_density = 0.5,
        pattern_fill = 'white',
        fill = 'lightgray',
        pattern_alpha = 0.3,
        alpha = 0.3
      )
  }
  
  if(facet_by_model) {
    cont_plot <- cont_plot + facet_grid(cols = vars(Model))
  }
  
  return(cont_plot)
}


# Data ####
## Diversity ####
data_diversity_plot <- get_div_plot(
  clones_over_time, 
  facet_by_model = FALSE,
  show_data = FALSE,
  add_pattern = FALSE,
  crop_y = FALSE
)

ggsave(
  filename = "img/clone_diversity.png",
  plot = data_diversity_plot,
  units = "cm",
  width = 10,
  height = 10
)

## SO ####
data_so_plot <- get_so_plot(
  so_occurence,
  facet_by_model = FALSE,
  show_data = FALSE,
  add_pattern = FALSE,
  crop_y = FALSE
)

ggsave(
  filename = "img/so_clones.png",
  plot = data_so_plot,
  units = "cm",
  width = 10,
  height = 10
)

## Contribution Z13264 ####
data_cont_Z13264_plot <- get_contribution_Z13264_plot(
  clone_contribution_dist,
  facet_by_model = FALSE,
  show_data = FALSE,
  add_pattern = FALSE,
  time_colors = c(
    "49" = "#E9727C",
    "300" = "#C21E2C",
    "1265" = "#7C131B"
  )
)

ggsave(
  filename = "img/clone_contribution_Z13264.png",
  plot = data_cont_Z13264_plot,
  units = "cm",
  width = 10,
  height = 10
)

## Contribution Z14004 ####
data_cont_Z14004_plot <- get_contribution_Z14004_plot(
  clone_contribution_dist,
  facet_by_model = FALSE,
  show_data = FALSE,
  add_pattern = FALSE,
  time_colors = c(
    "28" = "#24D2F9",
    "210" = "#04839F",
    "1111" = "#024150"
  )
)

ggsave(
  filename = "img/clone_contribution_Z14004.png",
  plot = data_cont_Z14004_plot,
  units = "cm",
  width = 10,
  height = 10
)

## Cell types ####
cell_types_plot <- ggplot(
  cell_by_group %>% filter(Cell_Group != "WBC"),
  aes(
    Time, Count,
    color = Cell_Group,
    group = Cell_Group
  )
) +
  geom_line(linewidth = 0.2) +
  geom_point(stroke = 0.16, size = 1) +
  geom_line(
    data = cell_by_group %>% filter(Cell_Group == "WBC"),
    # color = "black",
    linewidth = 0.4
  ) +
  geom_point(
    data = cell_by_group %>% filter(Cell_Group == "WBC"),
    # color = "black",
    shape = "cross",
    size = 1,
    stroke = 0.6
  ) +
  scale_y_log10() +
  scale_color_manual(values = c(
    brewer.pal(8, "Set2"),
    "black"
  )) +
  theme_classic() +
  labs(
    title = "Observed Cell Types",
    color = "Cell Type",
    x = "Time [days]",
    y = "Number of Cells"
  )

ggsave(
  filename = "img/cell_types.png",
  plot = cell_types_plot,
  units = "cm",
  width = 15,
  height = 8
)

# Model ####
## Generate Data ####
# get_best_sim_data <- function(ga_path, nr_sims = 1) {
#   source("src/optimization.R")
# 
#   load(ga_path)
# 
#   best_data <- run_one_cell(
#     ga_res@solution[1, ],
#     return_result = TRUE,
#     nr_sims = nr_sims,
#     sim_parallel = TRUE
#   )
# 
#   diversity_data <- bind_rows(
#     Z13264 = do.call(
#       rbind,
#       lapply(1:nr_sims, function(sim_idx)
#         get_diversity_df(best_data$Z13264[[sim_idx]]) %>%
#           mutate(Sim_Idx = sim_idx)
#       )
#     ),
#     Z14004 = do.call(
#       rbind,
#       lapply(1:nr_sims, function(sim_idx)
#         get_diversity_df(best_data$Z14004[[sim_idx]]) %>%
#           mutate(Sim_Idx = sim_idx)
#       )
#     ),
#     .id = "Subject"
#   )
# 
#   contribution_data <- bind_rows(
#     Z13264 = do.call(
#       rbind,
#       lapply(1:nr_sims, function(sim_idx)
#         get_contribution_dist_df(
#           best_data$Z13264[[sim_idx]],
#           break_list = break_list$Z13264
#         ) %>%
#           mutate(Sim_Idx = sim_idx)
#       )
#     ),
#     Z14004 = do.call(
#       rbind,
#       lapply(1:nr_sims, function(sim_idx)
#         get_contribution_dist_df(
#           best_data$Z14004[[sim_idx]],
#           break_list = break_list$Z14004
#         ) %>%
#           mutate(Sim_Idx = sim_idx)
#       )
#     ),
#     .id = "Subject"
#   )
# 
#   so_data <- bind_rows(
#     Z13264 = do.call(
#       rbind,
#       lapply(1:nr_sims, function(sim_idx)
#         get_so_df(best_data$Z13264[[sim_idx]]) %>%
#           mutate(Sim_Idx = sim_idx)
#       )
#     ),
#     Z14004 = do.call(
#       rbind,
#       lapply(1:nr_sims, function(sim_idx)
#         get_so_df(best_data$Z14004[[sim_idx]]) %>%
#           mutate(Sim_Idx = sim_idx)
#       )
#     ),
#     .id = "Subject"
#   )
# 
#   diff_data <- do.call(
#     rbind,
#     lapply(1:nr_sims, function(sim_idx) {
#       df <- get_full_diffs_df(list(
#         Z13264 = get_differences(list(best_data$Z13264[[sim_idx]]), 1, "Z13264"),
#         Z14004 = get_differences(list(best_data$Z14004[[sim_idx]]), 1, "Z14004")
#       )) %>%
#         mutate(Sim_Idx = sim_idx)
#       # %>%
#       #   merge(
#       #     point_weights,
#       #     all.x = TRUE,
#       #     all.y = FALSE
#       #   ) %>%
#       #   merge(
#       #     graph_weights,
#       #     by = c("Subject", "Graph"),
#       #     all.x = TRUE,
#       #     all.y = FALSE,
#       #     suffixes = c("_Point", "_Graph")
#       #   ) %>%
#       #   mutate(
#       #     Res =
#       #       Res*
#       #       coalesce(Weight_Point, 1)*
#       #       sqrt(coalesce(abs(Weight_Graph), 1))
#       #   ) %>%
#       #   mutate(Aspect = paste(Subject, Graph, sep = "_")) %>%
#       #   group_by(Aspect) %>%
#       #   mutate(RSS = sum(Res**2)) %>%
#       #   ungroup()
#       #
#       # df$Total <- sum((df %>% group_by(Aspect) %>% summarise(RSS = RSS[1]))$RSS)
#       #
#       return(df)
#     })
#   )
# 
#   return(list(
#     Diversity = diversity_data,
#     Contribution = contribution_data,
#     SO = so_data,
#     Diffs = diff_data
#   ))
# }

# best_two_compt <- get_best_sim_data(
#   "otm/data/two/final.Rda",
#   nr_sims = 100
# )
# 
# best_unified <- get_best_sim_data(
#   "otm/data/unified/final.Rda",
#   nr_sims = 100
# )
#
# best_one_compt <- get_best_sim_data(
#   "otm/data/one/final.Rda",
#   nr_sims = 100
# )
# 
# best_one_dw <- get_best_sim_data(
#   "otm/data/one_dw_SO/final.Rda",
#   nr_sims = 100
# )

# save(
#   best_two_compt, best_unified, best_one_compt, best_one_dw,
#   file = "data/paper_sim.Rda"
# )

load("data/paper_sim.Rda")

# One Compt. vs Two Compt. ####
## Diversity ####
div_plot <- get_div_plot(list(
    "One-compartment" = best_one_compt$Diversity,
    "Two-compartment" = best_two_compt$Diversity
  ))

ggsave(
  "img/clone_diversity_fit.png",
  plot = div_plot,
  units = "cm",
  width = 20,
  height = 10
)

## SO ####
so_plot <- get_so_plot(list(
  "One-compartment" = best_one_compt$SO,
  "Two-compartment" = best_two_compt$SO
))

ggsave(
  "img/so_clones_fit.png",
  plot = so_plot,
  units = "cm",
  width = 20,
  height = 10
)

so_one_zoom <- ggplot(
  best_one_compt$SO %>% filter(Time > 360),
  aes(
    Time, SO_Clone_Contribution,
    color = Subject,
    group = interaction(Subject)
  )
) +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun = "mean", geom = "point") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(
    values = color_scheme$Subject
  ) +
  labs(x = "", y = "")

ggsave(
  filename = "img/so_zoom.png",
  plot = so_one_zoom,
  units = "cm",
  width = 10,
  height = 6
)

## Contribution Z14004 ####
cont_Z14004_plot <- get_contribution_Z14004_plot(list(
  "One-compartment" = best_one_compt$Contribution,
  "Two-compartment" = best_two_compt$Contribution
))

ggsave(
  "img/clone_contribution_Z14004_fit.png",
  plot = cont_Z14004_plot,
  units = "cm",
  width = 20,
  height = 10
)

## Contribution Z13264 ####
cont_Z13264_plot <- get_contribution_Z13264_plot(list(
  "One-compartment" = best_one_compt$Contribution,
  "Two-compartment" = best_two_compt$Contribution
))

ggsave(
  "img/clone_contribution_Z13264_fit.png",
  plot = cont_Z13264_plot,
  units = "cm",
  width = 20,
  height = 10
)

## Residuals ####
### Adjust data ####
rename_graphs <- Vectorize(function(Graph_name) {
  if(grepl("Contribution_3", Graph_name)) {
    return("Clone Size Dist. #1")
  }
  
  if(grepl("Contribution_4", Graph_name)) {
    return("Clone Size Dist. #2")
  }
  
  if(grepl("SO", Graph_name)) {
    return("SO Clone Fraction")
  }
  
  return("Clone Diversity")
})

get_rs_df <- function(sim_data) {
  load("otm/weights.Rda")
  weight_df <- weight_df %>% rename(Subject = Animal_Id)
  
  return(merge(
    sim_data, 
    # point_weights,
    weight_df,
    all.x = TRUE,
    all.y = FALSE
  ) %>%
    mutate(
      Normalized_Res = Res*coalesce(Graph_Weight, 1)
    ) %>%
    group_by(Subject, Graph, Sim_Idx) %>%
    summarise(
      RSS = sum(Res**2),
      RSM = mean(Res**2),
      Normalized_RSS = sum(Normalized_Res**2),
      Normalized_RSM = mean(Normalized_Res**2),
      .groups = "drop"
    )
  )
}

rss_df <- bind_rows(
  "One-compartment" = get_rs_df(best_one_compt$Diffs),
  "Two-compartment" = get_rs_df(best_two_compt$Diffs),
  "Unified" = get_rs_df(best_unified$Diffs),
  "One-compartment (Downweighted SO)" = get_rs_df(best_one_dw$Diffs),
  .id = "Model"
)

mean_rss <- rss_df %>%
  filter(Model %in% c("One-compartment", "Two-compartment")) %>%
  group_by(Model, Graph, Subject) %>%
  summarise(
    Normalized_RSM = mean(Normalized_RSM),
    .groups = "drop"
  )

metric_order <- c(
  "SO",
  "Contribution_4",
  "Contribution_3",
  "Diversity"
)
progressive_df <- mean_rss %>%
  mutate(
    Graph = ifelse(grepl("Contribution_3", Graph), "Contribution_3", 
            ifelse(grepl("Contribution_4", Graph), "Contribution_4",
            Graph
            ))
  ) %>%
  pivot_wider(names_from = Graph, values_from = Normalized_RSM) %>%
  mutate(
    Contribution_3 = Diversity + Contribution_3,
    Contribution_4 = Contribution_3 + Contribution_4,
    SO = Contribution_4 + SO
  ) %>%
  pivot_longer(!c(Model, Subject), names_to = "Graph", values_to = "Normalized_RSM")
### Plot ####

rss_plot <- ggplot(
  mean_rss %>%
    mutate(Graph = rename_graphs(Graph), X = Graph),
  aes(
    interaction(sub("compartment", "comp.", Model), Subject, sep = "\n"),
    Normalized_RSM,
    group = Model,
    fill = Graph,
    pattern = Model
  )
) +
  geom_col_pattern(
    position = "dodge", color = "black",
    pattern_angle  = 45,
    pattern_spacing = 0.02
  ) +
  geom_col_pattern(
    color = "black", position = "dodge",
    pattern_angle  = 45,
    pattern_spacing = 0.02,
    data = progressive_df %>% 
      filter(Graph == metric_order[1]) %>%
      mutate(Graph = rename_graphs(Graph), X = "Total")
  ) +
  geom_col_pattern(
    color = "black", position = "dodge",
    pattern_angle  = 45,
    pattern_spacing = 0.02,
    data = progressive_df %>% 
      filter(Graph == metric_order[2]) %>%
      mutate(Graph = rename_graphs(Graph), X = "Total")
  ) +
  geom_col_pattern(
    color = "black", position = "dodge",
    pattern_angle  = 45,
    pattern_spacing = 0.02,
    data = progressive_df %>% 
      filter(Graph == metric_order[3]) %>%
      mutate(Graph = rename_graphs(Graph), X = "Total")
  ) +
  geom_col_pattern(
    color = "black", position = "dodge",
    pattern_angle  = 45,
    pattern_spacing = 0.02,
    data = progressive_df %>% 
      filter(Graph == metric_order[4]) %>%
      mutate(Graph = rename_graphs(Graph), X = "Total")
  ) +
  facet_wrap(vars(X), scales = "free") +
  theme_classic() +
  scale_fill_manual(values = color_scheme$Aspects) +
  theme(
    legend.position = "none"
  ) +
  labs(
    y = "Residuals",
    x = " "
  ) +
  scale_pattern_manual(values = c(
    "One-compartment" = "none",
    "Two-compartment" = "stripe"
  ))

ggsave(
  plot = rss_plot,
  filename = "img/rss.png",
  units = "cm",
  width = 25,
  height = 15
)

# Two Compt. vs Unified ####
## Diversity ####
div_plot <- get_div_plot(list(
  "Two-compartment" = best_two_compt$Diversity,
  "Unified" = best_unified$Diversity
))

ggsave(
  "img/clone_diversity_fit_unified.png",
  plot = div_plot,
  units = "cm",
  width = 20,
  height = 10
)

## SO ####
so_plot <- get_so_plot(list(
  "Two-compartment" = best_two_compt$SO,
  "Unified" = best_unified$SO
))

ggsave(
  "img/so_clones_fit_unified.png",
  plot = so_plot,
  units = "cm",
  width = 20,
  height = 10
)

## Contribution Z14004 ####
cont_Z14004_plot <- get_contribution_Z14004_plot(list(
  "Two-compartment" = best_two_compt$Contribution,
  "Unified" = best_unified$Contribution
))

ggsave(
  "img/clone_contribution_Z14004_fit_unified.png",
  plot = cont_Z14004_plot,
  units = "cm",
  width = 20,
  height = 10
)

## Contribution Z13264 ####
cont_Z13264_plot <- get_contribution_Z13264_plot(list(
  "Two-compartment" = best_two_compt$Contribution,
  "Unified" = best_unified$Contribution
))

ggsave(
  "img/clone_contribution_Z13264_fit_unified.png",
  plot = cont_Z13264_plot,
  units = "cm",
  width = 20,
  height = 10
)

# Downweighted One Compt. ####
## Diversity ####
div_plot <- get_div_plot(
  best_one_dw$Diversity, 
  facet_by_model = FALSE
)

ggsave(
  "img/clone_diversity_dw_one.png",
  plot = div_plot,
  units = "cm",
  width = 10,
  height = 10
)

## SO ####
so_plot <- get_so_plot(
  best_one_dw$SO, 
  facet_by_model = FALSE
)

ggsave(
  "img/so_clones_dw_one.png",
  plot = so_plot,
  units = "cm",
  width = 10,
  height = 10
)

## Contribution Z14004 ####
cont_Z14004_plot <- get_contribution_Z14004_plot(
  best_one_dw$Contribution, 
  facet_by_model = FALSE
)

ggsave(
  "img/clone_contribution_Z14004__dw_one.png",
  plot = cont_Z14004_plot,
  units = "cm",
  width = 10,
  height = 10
)

## Contribution Z13264 ####
cont_Z13264_plot <- get_contribution_Z13264_plot(
  best_one_dw$Contribution, 
  facet_by_model = FALSE
)

ggsave(
  "img/clone_contribution_Z13264_dw_one.png",
  plot = cont_Z13264_plot,
  units = "cm",
  width = 10,
  height = 10
)

# Rates ####
## Generate Data ####
load("otm/data/two/final.Rda")

best_rates <- ga_res@solution[1, c("pA", "dA", "tQA", "tAQ")]

best_fit_analitic <- derive_analitic_times(
  best_rates, kA = 1e5, kQ = 1e5
)

time_points <- c(
  10*c(1, 1.7, 3, 5.5),
  100*c(1, 1.7, 3, 5.5),
  1000*c(1, 1.7, 3, 5.5),
  10000*c(1, 1.7, 3, 5.5),
  100000
)

# empiric_measures <- measure_empiric_times(
#   100002,
#   1e5, 1e5, best_rates,
#   time_points = time_points
# )

load("data/empiric_time.Rda")

rename_action <- function(short_name) {
  switch(
    short_name,
    Act = "Activation",
    Deact = "Deactivation",
    Diff = "Differentiation",
    Prol = "Proliferation"
  )
}

empiric_time_summary <- bind_rows(lapply(rownames(empiric_measures$Counter), function(action_name) {
  qtl <- Hmisc::wtd.quantile(
    x = 1:max(time_points),
    weights = empiric_measures$Counter[action_name, 1:max(time_points)],
    probs = c(.25, .5, .75)
  )
  names(qtl) <- NULL
  mn <- weighted.mean(
    x = 1:max(time_points),
    w = empiric_measures$Counter[action_name, 1:max(time_points)]
  )
  
  data.frame(
    Action = rename_action(action_name),
    Stat = c("1st Qu.", "Median", "3rd Qu.", "Mean"),
    Value = c(qtl, mn)
  )
  
}))

empiric_prob_df <- bind_rows(lapply(rownames(empiric_measures$Counter), function(action_name) {
  final_time <- tail(which(empiric_measures$Counter[action_name, 1:max(time_points)] > 0), 1)
  total_cells <- sum(empiric_measures$Counter[action_name, 1:final_time])
  if(grepl("act", action_name, ignore.case = TRUE)) {
    
    space <- 200
    return(bind_rows(lapply(1:ceiling(final_time/space) - 1, function(curr_bin) {
      curr_tps <- ((curr_bin*space) + 1):((curr_bin+1)*space)
      curr_prob <- sum(empiric_measures$Counter[action_name, curr_tps])/total_cells
      
      return(data.frame(
        Action = rename_action(action_name),
        Time_Spaced = (curr_bin*space) + 1,
        Prob = curr_prob
      ))
    })) %>% filter(Prob > 0))
  } else {
    probs <- empiric_measures$Counter[action_name, 1:final_time]/total_cells
    names(probs) = NULL
    
    return(data.frame(
      Action = rename_action(action_name),
      Time_Spaced = 1:final_time,
      Prob = probs
    ) %>% filter(Prob > 0))
  }
}))

## Plots ####
best_fit_experimental <- empiric_measures$Measured_Times %>%
  t() %>% data.frame() %>% mutate(Time = time_points)
row.names(best_fit_experimental) <- NULL

avg_time_plot <- ggplot(
  best_fit_experimental %>%
    rename(
      Activation = Act,
      Deactivation = Deact,
      Proliferation = Prol,
      Differentiation = Diff
    ) %>%
    pivot_longer(!Time),
  aes(
    Time, value
  )
) +
  geom_point() +
  geom_hline(
    data = data.frame(t(best_fit_analitic)) %>%
      rename(
        Activation = Act,
        Deactivation = Deact,
        Proliferation = Prol,
        Differentiation = Diff
      ) %>%
      pivot_longer(c(
        Activation,
        Deactivation,
        Proliferation,
        Differentiation
      )),
    mapping = aes(yintercept = value),
    linewidth = 1,
    linetype = "dashed"
  ) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(vars(name), scales = "free_y") +
  theme_classic() +
  labs(
    x = "Simulation Time [days]",
    y = "Average Time to Event"
  )

ggsave(
  filename = "img/avg_time.png",
  plot = avg_time_plot,
  units = "cm",
  width = 20,
  height = 12
)

cells_per_time_plot <- ggplot(
  empiric_prob_df %>%
    filter(
      (Time_Spaced < 10000 & grepl("act", Action, ignore.case = TRUE)) |
      (Time_Spaced < 50 & !grepl("act", Action, ignore.case = TRUE))
    ),
  aes(Time_Spaced, Prob)
) +
  geom_vline(
    data = empiric_time_summary,
    mapping = aes(
      xintercept = Value,
      color = factor(Stat, levels = c("1st Qu.", "Median", "Mean", "3rd Qu."))
    ),
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point() +
  facet_wrap(vars(Action), scales = "free") +
  theme_classic() +
  labs(
    x = "Time [days]",
    y = "Probability",
    color = "Statistic"
  )

ggsave(
  filename = "img/cells_per_time.png",
  plot = cells_per_time_plot,
  units = "cm", 
  width = 20,
  heigh = 12
)

# AIC ####
## Adjust Data ####
get_n_pars <- function(model_name) {
  return(switch (model_name,
    "One-compartment" = 3,
    "One-compartment (Downweighted SO)" = 3,
    "Two-compartment" = 5,
    "Unified" = 6
  ))
}

n_obs <- (rss_df %>%
  filter(Model == "Unified" & Sim_Idx == 1) %>%
  mutate(n_obs = RSS/RSM))$n_obs %>%
  sum()

quality_metric_df <- rss_df %>%
  group_by(Model, Sim_Idx) %>%
  summarise(
    RSS = sum(Normalized_RSM),
    Sigma_Sq = sum(Normalized_RSM)/n_obs,
    N_Pars = get_n_pars(Model[1]),
    AIC = 2*N_Pars + n_obs*log(Sigma_Sq),
    BIC = log(n_obs)*N_Pars + n_obs*log(Sigma_Sq),
    .groups = "drop"
  )

## Plot ####
rss_dist_plt <- 
  ggplot(
  quality_metric_df %>% 
    filter(!grepl("Down", Model)) %>%
    pivot_longer(c(RSS)),
  aes(
    x = value, color = Model, y = after_stat(scaled)
  )
) +
  scale_x_log10() +
  geom_density(fill = "gray", alpha = 0.3, linewidth = 1) +
  # facet_grid(cols = vars(name), rows = vars(Model), scales = "free_x") +
  facet_wrap(vars(name), scales = "free_x") +
  labs(
    y = "Frequency",
    x = "Value"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text.y = element_text(size = 6)
  ) +
  scale_color_manual(values = c("#FF5555", "#33DD33", "#7777FF")) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

ggsave(
  filename = "img/rss_dist.png",
  plot = rss_dist_plt,
  units = "cm",
  width = 15, height = 10
)

aic_dist_plt <- ggplot(
  quality_metric_df %>% 
    filter(!grepl("Down", Model)) %>%
    pivot_longer(c(AIC)),
  aes(
    x = value, color = Model, y = after_stat(scaled)
  )
) +
  geom_density(fill = "gray", alpha = 0.3, linewidth = 1) +
  # facet_grid(cols = vars(name), rows = vars(Model), scales = "free_x") +
  facet_wrap(vars(name), scales = "free_x") +
  labs(
    y = "Frequency",
    x = "Value"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text.y = element_text(size = 6)
  ) +
  scale_color_manual(values = c("#FF5555", "#33DD33", "#7777FF")) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

ggsave(
  filename = "img/aic_dist.png",
  plot = aic_dist_plt,
  units = "cm",
  width = 15, height = 10
)

bic_dist_plt <- ggplot(
  quality_metric_df %>% 
    filter(!grepl("Down", Model)) %>%
    pivot_longer(c(BIC)),
  aes(
    x = value, color = Model, y = after_stat(scaled)
  )
) +
  geom_density(fill = "gray", alpha = 0.3, linewidth = 1) +
  # facet_grid(cols = vars(name), rows = vars(Model), scales = "free_x") +
  facet_wrap(vars(name), scales = "free_x") +
  labs(
    y = "Frequency",
    x = "Value"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.text.y = element_text(size = 6)
  ) +
  scale_color_manual(values = c("#FF5555", "#33DD33", "#7777FF")) +
  scale_y_continuous(breaks = c(0, 0.5, 1))

ggsave(
  filename = "img/bic_dist.png",
  plot = bic_dist_plt,
  units = "cm",
  width = 15, height = 12
)
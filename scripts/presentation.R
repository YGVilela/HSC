library(ggplot2)
library(tidyverse)

source("src/Simplified_HSC_Model.R")
source("src/utils.R")

load("data/extracted_radtke.Rda")
load("data/derived_radtke.Rda")

# Cell Types ####
get_cell_group <- function(cell_type) {
  group_names <- c(
    "Bcell",
    "Gran",
    "Lymph",
    "Mono",
    "NKcell",
    "Tcell",
    "WBC",
    "CD34+"
  )
  
  cell_groups <- rep("Other", length(cell_type))
  
  for(name in group_names) {
    cell_groups[which(grepl(name, cell_type))] <- name
  }
  
  return(cell_groups)
}

cell_by_group <- extracted_data %>%
  mutate(Cell_Group = get_cell_group(Cell_Type)) %>%
  group_by(Subject, Time, Cell_Group) %>%
  summarise(Count = sum(Count), .groups = "drop")

plt_cell_types <- ggplot(
  cell_by_group,
  aes(
    Time, Count,
    group = interaction(Subject, Cell_Group),
    fill = Cell_Group
  )
) + geom_col() +
  facet_grid(rows = vars(Subject)) +
  scale_y_log10() +
  labs(
    title = "Observed Cell Types",
    fill = "Cell Type"
  )

plt_cell_types
ggsave(
  filename = "img/cell_types.pdf",
  plot = plt_cell_types,
  device = "pdf",
  width = 15,
  height = 10,
  unit = "cm"
)

# Sample Counts ####
plt_observed_data <- ggplot(
  clone_data %>% 
    group_by(Subject, Time) %>%
    summarise(
      Cell_Count = sum(Clone_Size),
      Clone_Count = length(Clone_Size),
      .groups = "drop"
    ),
  aes(
    Time,
    group = Subject
  )
) + 
  geom_line(aes(y = Cell_Count, color = "# of Cells")) +
  geom_point(aes(y = Cell_Count, color = "# of Cells"), shape = "cross") +
  geom_line(aes(y = Clone_Count, color = "# of Clones")) +
  geom_point(aes(y = Clone_Count, color = "# of Clones"), shape = "cross") +
  facet_grid(rows = vars(Subject)) +
  scale_y_log10() + scale_x_log10() +
  labs(
    title = "Observed Cells",
    color = "Legend",
    y = "# of Cells/Clones"
  )

plt_observed_data
ggsave(
  filename = "img/plt_data.pdf",
  plot = plt_observed_data,
  device = "pdf",
  width = 15,
  height = 10,
  unit = "cm"
)

# Diversity (data) ####
plt_diversity <- ggplot(
  clones_over_time,
  aes(
    x = Time,
    y = Nr_Clones,
    color = Subject,
    group = Subject
  )
) +
  geom_line() +
  geom_point(shape = "cross") +
  labs(
    x = "Time",
    y = "Nr. Clones",
    title = "Clone Diversity",
    color = "Legend"
  )

plt_diversity
ggsave(
  filename = "img/plt_diversity.pdf", 
  plot = plt_diversity, 
  device = "pdf",
  width = 15,
  height = 10,
  units = "cm"
)

# Contribution (data) ####
plt_contribution <- ggplot(
  clone_contribution_dist %>% filter(
    Time %in% c(
      362, 397,
      462, 495
    )
  ),
  aes(
    x = Clone_Contribution,
    y = Frequency,
    group = interaction(Subject, Time),
    color = Subject
  )
) +
  geom_line() +
  geom_point(shape = "cross") +
  labs(
    x = "Clone Contribution",
    y = "Frequency",
    title = "Clone Contribution",
    color = "Legend"
  ) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(vars(paste0(Subject, "\nTime = ", Time)))

plt_contribution
ggsave(
  filename = "img/plt_contribution.pdf", 
  plot = plt_contribution, 
  device = "pdf",
  width = 25,
  height = 15,
  units = "cm"
))()

# SO Clones ####
plt_so <- ggplot(
  so_occurence %>%
    filter(Total_Cells > 1000),
  aes(
    x = Time,
    y = SO_Clone_Contribution,
    color = Subject,
    group = Subject
  )
) +
  geom_line() +
  geom_point(shape = "cross") +
  labs(
    x = "Time",
    y = "SO Clone Contribution",
    title = "Single Occurring Clones",
    color = "Legend"
  )

plt_so
ggsave(
  filename = "img/plt_so.pdf", 
  plot = plt_so, 
  device = "pdf",
  width = 15,
  height = 10,
  units = "cm"
)

# Diversity (One Compartment Fit) ####
nr_clones <- 17e4
pA = 0.7535011
dA = 0.1
one_comp_sim <- simulate_HSC(
  nr_clones = nr_clones,
  pA = pA, dA = dA,
  tQA = 0, tAQ = 0,
  pQ = 0,
  kA = nr_clones*pA/(pA - dA), kQ = 1,
  init_A = 1, init_Q = 0,
  tps = all_data_tps
)

one_comp_div <- get_diversity_df(one_comp_sim)
plt_one_comp_div <- ggplot(
  clones_over_time,
  aes(
    x = Time,
    y = Nr_Clones,
    color = Subject,
    group = Subject
  )
) +
  geom_line(
    data = one_comp_div,
    mapping = aes(
      color = NULL,
      group = NULL
    ),
    linetype = "solid"
  ) +
  geom_line(linetype = "dashed") +
  geom_point(shape = "cross") +
  labs(
    x = "Time",
    y = "Nr. Clones",
    title = "Clone Diversity Fit",
    color = "Legend"
  ) +
  ylim(0, 35000)

plt_one_comp_div
ggsave(
  filename = "img/plt_one_comp_div.pdf", 
  plot = plt_one_comp_div, 
  device = "pdf",
  width = 15,
  height = 10,
  units = "cm"
)

# Contribution (One compartment fit) ####
one_comp_cont <- bind_rows(
  Z13264 = get_contribution_dist_df(one_comp_sim[, c("t397", "t462")]),
  Z14004 = get_contribution_dist_df(one_comp_sim[, c("t362", "t495")]),
  .id = "Subject"
)

plt_one_comp_cont <- ggplot(
  clone_contribution_dist %>% filter(
    Time %in% c(
      362, 397,
      462, 495
    )
  ),
  aes(
    x = Clone_Contribution,
    y = Frequency,
    group = interaction(Subject, Time),
    color = Subject
  )
) +
  geom_line(
    data = one_comp_cont,
    mapping = aes(color = NULL)
  ) +
  geom_line(linetype = "dashed") +
  geom_point(shape = "cross") +
  labs(
    x = "Clone Contribution",
    y = "Frequency",
    title = "Clone Contribution",
    color = "Legend"
  ) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(vars(paste0(Subject, "\nTime = ", Time)))

plt_one_comp_cont
ggsave(
  filename = "img/plt_one_comp_cont.pdf", 
  plot = plt_one_comp_cont, 
  device = "pdf",
  width = 25,
  height = 15,
  units = "cm"
)

# Sample Effect on contribution ####
sample_effect_df <- get_contribution_dist_df(cbind(
  "t10" = sample_population(one_comp_sim[, "t397"], floor(sum(one_comp_sim[, "t397"])*0.1)),
  "t25" = sample_population(one_comp_sim[, "t397"], floor(sum(one_comp_sim[, "t397"])*0.25)),
  "t100" = one_comp_sim[, "t397"]
))

plt_sample_effect <- ggplot(
  sample_effect_df,
  aes(
    x = Clone_Contribution,
    y = Frequency,
    group = Time
  )
) +
  geom_line() +
  scale_y_log10() + scale_x_log10() +
  facet_grid(cols = vars(
    factor(
      paste0("Sampling ", Time, "%"),
      levels = c(
        "Sampling 100%", "Sampling 25%", "Sampling 10%"
      )
    )
  )) +
  labs(
    title = "Sampling Effect",
    x = "Clone Contribution"
  )

plt_sample_effect
ggsave(
  filename = "img/plt_sample_effect.pdf", 
  plot = plt_sample_effect, 
  device = "pdf",
  width = 25,
  height = 10,
  units = "cm"
)

# Contribution (One compartment fit with sample) ####
nr_clones <- floor(16.07755e4)
pA = 0.7535011
dA = 0.1328183
one_comp_sample_sim_Z13264 <- simulate_HSC(
  nr_clones = nr_clones,
  pA = pA, dA = dA,
  tQA = 0, tAQ = 0,
  pQ = 0,
  kA = 1.083431*nr_clones*pA/(pA - dA), kQ = 1,
  init_A = 1, init_Q = 0,
  tps = data_tps$Z13264,
  sample_sizes = sample_sizes_vector$Z13264
)

one_comp_sample_sim_Z14004 <- simulate_HSC(
  nr_clones = nr_clones,
  pA = pA, dA = dA,
  tQA = 0, tAQ = 0,
  pQ = 0,
  kA = nr_clones*pA/(pA - dA), kQ = 1,
  init_A = 1, init_Q = 0,
  tps = data_tps$Z14004,
  sample_sizes = sample_sizes_vector$Z14004
)

one_comp_sample_cont <- bind_rows(
  Z13264 = get_contribution_dist_df(one_comp_sample_sim_Z13264, break_list = break_list$Z13264),
  Z14004 = get_contribution_dist_df(one_comp_sample_sim_Z14004, break_list = break_list$Z14004),
  .id = "Subject"
)

plt_one_comp_sample_cont <- ggplot(
  clone_contribution_dist %>% filter(
    Time %in% c(
      362, 397,
      462, 495
    )
  ),
  aes(
    x = Clone_Contribution,
    y = Frequency,
    group = interaction(Subject, Time),
    color = Subject
  )
) +
  geom_line(
    data = one_comp_sample_cont %>% filter(
      Time %in% c(
        362, 397,
        462, 495
      )
    ),
    mapping = aes(color = NULL)
  ) +
  geom_line(linetype = "dashed") +
  geom_point(shape = "cross") +
  labs(
    x = "Clone Contribution",
    y = "Frequency",
    title = "Clone Contribution",
    color = "Legend"
  ) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(vars(paste0(Subject, "\nTime = ", Time)))

plt_one_comp_sample_cont
ggsave(
  filename = "img/plt_one_comp_sample_cont.pdf", 
  plot = plt_one_comp_sample_cont, 
  device = "pdf",
  width = 25,
  height = 15,
  units = "cm"
)

one_comp_sample_div <- bind_rows(
  Z13264 = get_diversity_df(one_comp_sample_sim_Z13264),
  Z14004 = get_diversity_df(one_comp_sample_sim_Z14004),
  .id = "Subject"
)

plt_one_comp_sample_div <- ggplot(
  clones_over_time,
  aes(
    x = Time,
    y = Nr_Clones,
    color = Subject,
    group = Subject
  )
) +
  geom_line(
    data = one_comp_sample_div,
    aes(color = NULL),
    linetype = "solid",
    color = "black"
  ) +
  geom_line(linetype = "dashed") +
  geom_point(shape = "cross") +
  labs(
    x = "Time",
    y = "Nr. Clones",
    title = "Clone Diversity Fit",
    color = "Legend"
  ) +
  facet_wrap(vars(Subject))

plt_one_comp_sample_div
ggsave(
  filename = "img/plt_one_comp_sample_div.pdf",
  plot = plt_one_comp_sample_div, 
  device = "pdf",
  width = 20,
  height = 10,
  units = "cm"
)

one_comp_sample_so <- bind_rows(
  Z13264 = get_so_df(one_comp_sample_sim_Z13264),
  Z14004 = get_so_df(one_comp_sample_sim_Z14004),
  .id = "Subject"
)

plt_one_comp_sample_so <- ggplot(
  so_occurence %>%
    filter(Total_Cells > 1000),
  aes(
    x = Time,
    y = SO_Clone_Contribution,
    color = Subject,
    group = Subject
  )
) +
  geom_line(
    data = one_comp_sample_so,
    mapping = aes(
      y = SO_Clone_Contribution,
      color = NULL,
      group = NULL
    )
  ) +
  geom_line(linetype = "dashed") +
  geom_point(shape = "cross") +
  labs(
    x = "Time",
    y = "SO Clone Contribution",
    title = "Single Occurring Clones",
    color = "Legend"
  ) +
  facet_wrap(vars(Subject))

plt_one_comp_sample_so
ggsave(
  filename = "img/plt_one_comp_sample_so.pdf",
  plot = plt_one_comp_sample_so, 
  device = "pdf",
  width = 20,
  height = 10,
  units = "cm"
)

# Two compartment fits ####
load("data/one_cell/ga_two_final.Rda")

best_ga <- run_one_cell(ga_res@solution[1, ], return_result = TRUE)

two_comp_so <- bind_rows(
  Z13264 = get_so_df(best_ga$Z13264[[1]]),
  Z14004 = get_so_df(best_ga$Z14004[[1]]),
  .id = "Subject"
)

plt_two_comp_so <- ggplot(
  so_occurence %>%
    filter(Total_Cells > 1000),
  aes(
    x = Time,
    y = SO_Clone_Contribution,
    color = Subject,
    group = Subject
  )
) +
  geom_line(
    data = two_comp_so %>%
      filter(Total_Cells > 1000),
    mapping = aes(
      y = SO_Clone_Contribution,
      color = NULL,
      group = NULL
    )
  ) +
  geom_line(linetype = "dashed") +
  geom_point(shape = "cross") +
  labs(
    x = "Time",
    y = "SO Clone Contribution",
    title = "Single Occurring Clones",
    color = "Legend"
  ) +
  facet_wrap(vars(Subject))

plt_two_comp_so
ggsave(
  filename = "img/plt_two_comp_so.pdf",
  plot = plt_two_comp_so, 
  device = "pdf",
  width = 20,
  height = 10,
  units = "cm"
)

two_comp_div <- bind_rows(
  Z13264 = get_diversity_df(best_ga$Z13264[[1]]),
  Z14004 = get_diversity_df(best_ga$Z14004[[1]]),
  .id = "Subject"
)

plt_two_comp_div <- ggplot(
  clones_over_time,
  aes(
    x = Time,
    y = Nr_Clones,
    color = Subject,
    group = Subject
  )
) +
  geom_line(
    data = two_comp_div,
    aes(color = NULL),
    linetype = "solid",
    color = "black"
  ) +
  geom_line(linetype = "dashed") +
  geom_point(shape = "cross") +
  labs(
    x = "Time",
    y = "Nr. Clones",
    title = "Clone Diversity Fit",
    color = "Legend"
  ) +
  facet_wrap(vars(Subject))

plt_two_comp_div
ggsave(
  filename = "img/plt_two_comp_div.pdf",
  plot = plt_two_comp_div, 
  device = "pdf",
  width = 20,
  height = 10,
  units = "cm"
)

two_comp_cont <- bind_rows(
  Z13264 = get_contribution_dist_df(best_ga$Z13264[[1]]),
  Z14004 = get_contribution_dist_df(best_ga$Z14004[[1]]),
  .id = "Subject"
)

plt_two_comp_cont <- ggplot(
  clone_contribution_dist %>% filter(
    Time %in% c(
      362, 397,
      462, 495
    )
  ),
  aes(
    x = Clone_Contribution,
    y = Frequency,
    group = interaction(Subject, Time),
    color = Subject
  )
) +
  geom_line(
    data = two_comp_cont %>% filter(
      Time %in% c(
        362, 397,
        462, 495
      )
    ),
    mapping = aes(color = NULL)
  ) +
  geom_line(linetype = "dashed") +
  geom_point(shape = "cross") +
  labs(
    x = "Clone Contribution",
    y = "Frequency",
    title = "Clone Contribution",
    color = "Legend"
  ) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(vars(paste0(Subject, "\nTime = ", Time)))

plt_two_comp_cont
ggsave(
  filename = "img/plt_two_comp_cont.pdf", 
  plot = plt_two_comp_cont, 
  device = "pdf",
  width = 25,
  height = 15,
  units = "cm"
)

diffs_df <- list(
  Z13264 = get_differences(best_ga$Z13264, 1, "Z13264"),
  Z14004 = get_differences(best_ga$Z14004, 1, "Z14004")
) %>% get_full_diffs_df()

normalized_res <-
  merge(
    diffs_df, point_weights,
    all.x = TRUE,
    all.y = FALSE
  ) %>%
  mutate(
    Normalized_Res = Res*coalesce(Weight, 1)
  )

ggplot(normalized_res, aes(Normalized_Res, after_stat(density))) +
  geom_histogram(color = "blue") +
  geom_density(color = "blue")

ggplot(normalized_res, aes(x = Normalized_Res, after_stat(density), group = interaction(Subject, Graph), color = paste0(Subject, "; ", Graph))) +
  geom_histogram() +
  geom_density() +
  facet_grid(rows = vars(paste0(Subject, "\n", Graph)), scales = "free")
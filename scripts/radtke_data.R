library(readr)
library(tidyverse)
source("src/utils.R")

# Read tables ####
# Tables available on https://github.com/KiemLab-RIS/ISA-Clone-tracking
Z13264_globalAlignments <- readr::read_table(
  file = "data/Z13264_globalAlignments.txt.gz",
  col_names = c(
    "Count",
    "Bar_Code",
    "Description",
    "Var1",
    "Sequence"
  ), 
  col_types = c(
    col_double(),
    col_character(),
    col_character(),
    col_double(),
    col_character()
  )
)

Z14004_globalAlignments <- readr::read_table(
  file = "data/Z14004_globalAlignments.txt.gz",
  col_names = c(
    "Count",
    "Bar_Code",
    "Description",
    "Var1",
    "Sequence"
  ), 
  col_types = c(
    col_double(),
    col_character(),
    col_character(),
    col_double(),
    col_character()
  )
)

raw_data <- rbind(
  Z13264_globalAlignments,
  Z14004_globalAlignments
)

# Extract data ####
clone_id <- raw_data[
  is.na(raw_data$Description),
  c("Count", "Bar_Code")
]
clone_id <- separate(
  data = clone_id,
  col = "Count",
  into = c("Subject", "Clone_Id"), sep = "_",
  extra = "merge",
  remove = TRUE
)
clone_id$Bar_Code <- substring(clone_id$Bar_Code, 2)

read_data <- raw_data[
  !is.na(raw_data$Description),
  c("Count", "Bar_Code", "Description")
]
read_data <- separate(
  data = read_data,
  col = "Description",
  into = c("Subject", "Time", "Source", "Cell_Type"),
  sep = "_",
  extra = "merge",
  remove = FALSE
)
read_data$Count <- as.numeric(read_data$Count)
read_data$Time <- as.numeric(sub("DPT", "", read_data$Time))
read_data$Cell_Type <- sub(".tran", "", read_data$Cell_Type)
#read_data <- drop_na(read_data)

extracted_data <- merge(clone_id, read_data, all = TRUE)
extracted_data <- drop_na(extracted_data)

# Clone data ####
clone_data <- extracted_data %>%
  filter(Source == "PB" & grepl("WBC", Cell_Type)) %>%
  group_by(Subject, Time, Bar_Code) %>%
  summarize(Clone_Size = sum(Count), .groups = "drop_last") %>%
  mutate(
    Cell_Count = sum(Clone_Size),
    Mean_Clone_Size = mean(Clone_Size),
    Clone_Contribution = Clone_Size/sum(Clone_Size)
  ) %>%
  ungroup()

# Clone distribution ####

clone_contribution_dist <- clone_data %>%
  group_by(Subject, Time) %>%
  group_modify(function(curr_data, keys) {
    hist_data <- hist(log10(curr_data$Clone_Size/curr_data$Cell_Count), plot = FALSE)
    non_zero <- which(hist_data$counts > 0)
    
    return(data.frame(
      Clone_Contribution = 10**hist_data$breaks[non_zero],
      Frequency = hist_data$counts[non_zero]/length(curr_data$Clone_Size)
    ))
  }) %>%
  ungroup()

clone_contribution_dist$Engrafment_Phase <- factor(
  engrafment_phase(clone_contribution_dist$Time),
  levels = c(
    "Initial Engrafment",
    "Stabilization",
    "Homeostasis",
    "Long-term Engrafment"
  )
)

# Clones over Time ####
clones_over_time <- clone_data %>%
  group_by(Subject) %>%
  group_modify(function(curr_data, keys) {
    tps <- unique(curr_data$Time)
    
    nr_clones <- sapply(tps, function(tp) 
      curr_data %>%
        filter(Time >= tp) %>%
        select(Bar_Code) %>%
        unique() %>% nrow()
    )
    
    return(data.frame(
      Time = tps,
      Nr_Clones = nr_clones
    ))
  }) %>%
  ungroup()

# Single Occurring Clones ####
so_occurence <- extracted_data %>%
  filter(regexpr("WBC", Cell_Type) > 0 & Source == "PB") %>%
  group_by(Subject, Bar_Code) %>%
  mutate(Occurrences = length(unique(Time))) %>%
  group_by(Subject, Time) %>%
  group_modify(function(curr_data, keys) {
    so_clones <- length(unique((curr_data %>% filter(Occurrences == 1))$Bar_Code))
    total_clones = length(unique(curr_data$Bar_Code))
    so_cells = sum((curr_data %>% filter(Occurrences == 1))$Count)
    total_cells = sum(curr_data$Count)
    
    return(data.frame(
      SO_Cells = so_cells,
      Total_Cells = total_cells,
      SO_Cell_Contribution = so_cells/total_cells,
      
      SO_Clones = so_clones,
      Total_Clones = total_clones,
      SO_Clone_Contribution = so_clones/total_clones
    ))
  }) %>%
  ungroup()


# Other data ####
data_tps <- list(
  "Z13264" = (clones_over_time %>% filter(Subject == "Z13264"))$Time,
  "Z14004" = (clones_over_time %>% filter(Subject == "Z14004"))$Time
)
all_data_tps <- sort(clones_over_time$Time)

data_sample_sizes <- list(
  "Z13264" = clone_data %>%
    filter(Subject == "Z13264") %>%
    group_by(Time) %>%
    summarize(Cell_Count = min(Cell_Count), .groups = "drop"),
  
  "Z14004" = clone_data %>%
    filter(Subject == "Z14004") %>%
    group_by(Time) %>%
    summarize(Cell_Count = min(Cell_Count), .groups = "drop")
)

sample_sizes_vector <- list(
  "Z13264" = data_sample_sizes[["Z13264"]]$Cell_Count,
  "Z14004" = data_sample_sizes[["Z14004"]]$Cell_Count
)
names(sample_sizes_vector[["Z13264"]]) <- paste0("t", data_tps[["Z13264"]])
names(sample_sizes_vector[["Z14004"]]) <- paste0("t", data_tps[["Z14004"]])

# Get contributions histogram breaks
break_list <- list(
  "Z13264" = lapply(data_tps[["Z13264"]], function(tp) {
    curr_data <- clone_data %>%
      filter(Subject == "Z13264" & Time == tp)
    
    hist_data <- hist(log10(curr_data$Clone_Size/curr_data$Cell_Count), plot = FALSE)
    
    return(c(-Inf, hist_data$breaks, 0))
  }),
  
  "Z14004" = lapply(data_tps[["Z14004"]], function(tp) {
    curr_data <- clone_data %>%
      filter(Subject == "Z14004" & Time == tp)
    
    hist_data <- hist(log10(curr_data$Clone_Size/curr_data$Cell_Count), plot = FALSE)
    
    return(c(-Inf, hist_data$breaks, 0))
  })
)

names(break_list[["Z13264"]]) <- paste0("t", data_tps[["Z13264"]])
names(break_list[["Z14004"]]) <- paste0("t", data_tps[["Z14004"]])

# Save data ####
save(
  raw_data, extracted_data,
  file = "data/extracted_radtke.Rda"
)
save(
  clone_data,
  clone_contribution_dist,
  clones_over_time,
  so_occurence,
  
  data_tps,
  all_data_tps,
  sample_sizes_vector,
  break_list,
  
  file = "data/derived_radtke.Rda"
)

# Plots ####
ggplot(
  data = clones_over_time,
  aes(
    Time, Nr_Clones,
    color = Subject,
    group = Subject
  )
) +
  geom_point(shape = "cross", stroke = 1) +
  geom_line() +
  labs(
    y = "Nr. Clones",
    title = "Clone Diversity"
  )

ggsave(
  filename = "img/clone_diversity.png",
  units = "cm",
  width = 15,
  height = 10
)

ggplot(
  data = so_occurence %>%
    filter(Total_Cells > 1000),
  aes(
    Time, SO_Clone_Contribution,
    color = Subject,
    group = Subject
  )
) +
  geom_point(shape = "cross", stroke = 1) +
  geom_line() +
  labs(
    y = "Fraction of SO Clones",
    title = "SO Clones"
  )

ggsave(
  filename = "img/so_occurrence.png",
  units = "cm",
  width = 15,
  height = 10
)

ggplot(
  data = clone_contribution_dist %>%
    filter(Time %in% c(362,397,462,495)),
  aes(
    Clone_Contribution, Frequency,
    color = Subject,
    group = Subject
  )
) +
  geom_point(shape = "cross", stroke = 1) +
  geom_line() +
  facet_wrap(vars(Time)) +
  scale_y_log10() + scale_x_log10() +
  labs(
    y = "Fraction of SO Clones",
    title = "Clone Contribution"
  )

ggsave(
  filename = "img/clone_contribution.png",
  units = "cm",
  width = 25,
  height = 15
)

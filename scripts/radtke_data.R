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
      Frequency = hist_data$counts[non_zero]/curr_data$Cell_Count[1]
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

# Infer clone size data based on future reads ####
inferred_clone_data <- clone_data %>%
  group_by(Subject) %>%
  group_map(function(curr_data, keys) {
    tps <- sort(unique(curr_data$Time))
    clones <- unique(curr_data$Bar_Code)
    
    infered_data <- matrix(
      data = 0,
      nrow = length(clones),
      ncol = length(tps) + 1,
      dimnames = list(
        clones,
        paste0("t", c(0, tps))
      )
    )
    infered_data[, "t0"] <- 1
    
    pb <- progress::progress_bar$new(
      format = "[:bar :current/:total] (:eta)",
      total = length(clones)
    )
    for(curr_clone in clones) {
      curr_clone_data <- curr_data %>% filter(Bar_Code == curr_clone)
      
      # Get clone appearances over time
      appearances <- sort(curr_clone_data$Time)
      
      last_appearance_idx <- 0
      last_appearance <- 0
      last_size <- 1
      
      final_appearence <- max(appearances)
      
      # For each time point until its final appearance
      for(tp in tps[which(tps <= final_appearence)]) {
        
        # If it appears on that read, use that data
        if(tp %in% appearances) {
          last_size <- (curr_clone_data %>% filter(Time == tp))$Clone_Size
          last_appearance <- tp
          last_appearance_idx <- last_appearance_idx + 1
          
          infered_data[curr_clone, paste0("t", tp)] <- last_size
        }
        # Otherwise, complete based on last and next appearance assuming linear growth
        # TODO: If we decide to fill the data, consider using another growth (?)
        else {
          next_appearance <- appearances[last_appearance_idx + 1]
          next_size <- (curr_clone_data %>% filter(Time == next_appearance))$Clone_Size
          infered_data[curr_clone, paste0("t", tp)] <- floor(
            (last_size*(next_appearance - tp) + next_size*(tp - last_appearance))/(next_appearance - last_appearance)
          )
        }
      }
      pb$tick()
    }
    
    return(infered_data)
  })

names(inferred_clone_data) <- c("Z13264", "Z14004")

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
  inferred_clone_data,
  file = "data/derived_radtke.Rda"
)

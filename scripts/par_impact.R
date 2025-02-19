source("src/optimization.R")

load("data/one_cell_v2/ga_two_final.Rda")
base_pars <- ga_res@solution[1,]

get_graph_data <- function(curr_pars) {
  
  res <- simulate_HSC(
    nr_clones = curr_pars["nr_clones"],
    pA = curr_pars["pA"],
    dA = curr_pars["dA"],
    tQA = curr_pars["tQA"],
    tAQ = curr_pars["tAQ"],
    pQ = 0,
    kA = curr_pars["kA"],
    kQ = curr_pars["kQ"],
    init_A = 1, init_Q = 0,
    tps = data_tps$Z14004,
    sample_sizes = sample_sizes_vector$Z14004
  )
  
  div <- get_diversity_df(res) %>%
    mutate(
      X = Time, 
      Y = Nr_Clones,
      Aspect = "Diversity",
      .keep = "none"
    )
  
  so <- get_so_df(res) %>%
    mutate(
      X = Time, 
      Y = SO_Clone_Contribution,
      Aspect = "SO",
      .keep = "none"
    )
  
  cont <- get_contribution_dist_df(
      res, 
      break_list = break_list$Z14004
    ) %>%
      mutate(
        X = Clone_Contribution,
        Y = Frequency,
        Aspect = paste("Contribution", Time, sep = "_"),
        .keep = "none"
      )
  
  return(bind_rows(div, so, cont))
}

evaluated_pars <- list(
  pA = runif(1, 0.2, 0.9),
  dA = runif(1, 0.025, 0.15),
  tQA = runif(1, 0.0002, 0.02),
  tAQ = runif(1, 0.002, 0.02),
  nr_clones = runif(1, 1e4, 1e5),
  kA = runif(1, 6e4, 6e5),
  kQ = runif(1, 6e4, 6e5)
)

evaluated_data <- do.call(rbind, lapply(names(evaluated_pars), function(par_name) {
  par_list <- evaluated_pars[[par_name]]
  
  sim_data <- do.call(rbind, parallel::mclapply(par_list, function(par_value) {
    load("data/one_cell_v2/ga_pQ_final.Rda")
    base_pars <- ga_res@solution[1,]
    
    base_pars["nr_clones"] <- floor(base_pars["clone_mult"]*1e4)
    base_pars["kA"] <- with(as.list(base_pars), {
      eq_A <- space_mult_A*nr_clones
      
      return(
        eq_A*pA/(pA - dA)
      )
    })
    base_pars["kQ"] <- with(as.list(base_pars), {
      eq_A <- space_mult_A*nr_clones
      eq_Q <- space_mult_Q*nr_clones
      
      return(ifelse(
        tAQ > 0 & tQA > 0,
        eq_Q*(tAQ*pA*eq_A)/(tAQ*pA*eq_A - tQA*dA*eq_Q),
        1
      ))
    })
    
    curr_pars <- base_pars
    curr_pars[par_name] <- par_value
    
    print(paste("Starting", par_name, "=", par_value))
    st <- Sys.time()
    res <-
      get_graph_data(curr_pars) %>%
      mutate(
        Par_Name = par_name,
        Par_Value = par_value
      )
    et <- Sys.time()
    print(paste("Done with", par_name, "=", par_value))
    print(et - st)  
    return(res)
  }, mc.cores = min(length(par_list), floor(parallel::detectCores()/2) - 1)))
  
  return(sim_data)
}))
evaluated_data$X <- as.numeric(evaluated_data$X)
evaluated_data$Y <- as.numeric(evaluated_data$Y)
evaluated_data$Par_Value <- as.numeric(evaluated_data$Par_Value)
evaluated_data <- drop_na(evaluated_data)

save(
  evaluated_data,
  file = "data/par_impact.Rda"
)

stop("Done")
library(ggplot2)
library(tidyr)
library(dplyr)

load("data/par_impact.Rda")

# Diversity
plt_div <- ggplot(
  evaluated_data %>%
    filter(
      Aspect == "Diversity"
    ) %>%
    mutate(
      Stage = engrafment_phase(X)
    ) %>%
    filter(Stage %in% c("Long-term Engrafment")) %>%
    group_by(Par_Name, Par_Value, Stage) %>%
    summarise(
      Y = mean(Y),
      .groups = "drop"
    ),
  aes(
    x = Par_Value,
    y = Y,
    group = interaction(Par_Name, Stage),
    # color = Stage
  )
) +
  geom_point() +
  geom_smooth() +
  facet_wrap(vars(Par_Name), scales = "free_x") +
  labs(
    title = "Parameter Impact on Late Diversity",
    x = "Par. Value",
    y = "Nr. Clones"
  )

ggsave(
  plot = plt_div,
  filename = "img/par_impact_diversity.png",
  units = "cm",
  width = 20,
  height = 20
)

# SO
plt_so <- ggplot(
  evaluated_data %>%
    filter(
      Aspect == "SO"
    ) %>%
    mutate(
      Stage = engrafment_phase(X)
    ) %>%
    filter(Stage %in% c("Long-term Engrafment")) %>%
    group_by(Par_Name, Par_Value, Stage) %>%
    mutate(
      Y = mean(Y)
    ),
  aes(
    x = Par_Value,
    y = Y,
    group = interaction(Par_Name, Stage)
  )
) +
  geom_point() +
  geom_smooth() +
  facet_wrap(vars(Par_Name), scales = "free_x") +
  labs(
    title = "Parameter Impact on Late SO Clones",
    x = "Par. Value",
    y = "Fraction of SO Clones"
  )

ggsave(
  plot = plt_so,
  filename = "img/par_impact_so.png",
  units = "cm",
  width = 20,
  height = 20
)

# Contribution
plt_cont <- ggplot(
  evaluated_data %>%
    filter(
      Aspect == "Contribution_495"
    ) %>%
    mutate(
      Stage = ifelse(X < 10**(-4.5), "Small", ifelse(X > 10**(-3.8) & X < 1e-3, "Big", ifelse(X >= 1e-3, "Huge", "Mid")))
    ) %>% 
    filter(Stage %in% c("Small", "Big")) %>%
    group_by(Par_Name, Par_Value, Stage) %>%
    mutate(
      Y = mean(Y)
    ),
  aes(
    x = Par_Value,
    y = Y,
    group = interaction(Par_Name, Stage),
  )
) +
  geom_point(aes(shape = Stage), stroke = 1) +
  geom_smooth(aes(color = Stage)) +
  facet_wrap(vars(Par_Name), scales = "free_x") +
  labs(
    title = "Parameter Impact on Clone Contribution",
    x = "Par. Value",
    y = "Frequency"
  ) +
  scale_color_manual(
    values = c(
      "Small" = "blue",
      "Big" = "green"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Small" = "bullet",
      "Big" = "cross"
    )
  )

ggsave(
  plot = plt_cont,
  filename = "img/par_impact_contribution.png",
  units = "cm",
  width = 20,
  height = 20
)
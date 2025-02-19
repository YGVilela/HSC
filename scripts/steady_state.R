source("src/optimization.R")
load("data/one_cell/ga_pQ_final.Rda")
test_pQ_1 <- get_one_cell_actions(
  ga_res@solution[1, ],
  nr_sims = 1,
  sim_parallel = TRUE
)

load("data/one_cell_v2/ga_pQ_final.Rda")
test_pQ_2 <- get_one_cell_actions(
  ga_res@solution[1, ],
  nr_sims = 1,
  sim_parallel = TRUE
)

load("data/one_cell/ga_two_final.Rda")
test_two_1 <- get_one_cell_actions(
  ga_res@solution[1, ],
  nr_sims = 1,
  sim_parallel = TRUE
)

load("data/one_cell_v2/ga_two_final.Rda")
test_two_2 <- get_one_cell_actions(
  ga_res@solution[1, ],
  nr_sims = 1,
  sim_parallel = TRUE
)

with(as.list(ga_res@solution[1, ]), {
  nr_clones <- ceiling(clone_mult*1e4)

  eq_A <- space_mult_A*nr_clones
  eq_Q <- space_mult_Q*nr_clones
  
  kA <- eq_A*pA/(pA - dA)
  
  kQ <- eq_Q*(tAQ*pA*eq_A)/(tAQ*pA*eq_A - tQA*dA*eq_Q)
  
  alpha <- (pA - dA)/pA
  
  # res <- 1/tAQ + (kA/kQ)*(pA/dA - 1)*(1/tQA)
  res <- 1/tAQ + (kA/kQ)*(pA/(dA*tQA) - 1/tQA)
    
  print(res)
})

bind_rows(
  # "With pQ (06.01.25)" = colMeans(t(sapply(1:1, function(idx) rowMeans(test_pQ_1[[idx]][, 250:1000])))),
  # "With pQ (12.01.25)" = colMeans(t(sapply(1:1, function(idx) rowMeans(test_pQ_2[[idx]][, 250:1000])))),
  "No pQ (06.01.25)" = colMeans(t(sapply(1:1, function(idx) rowMeans(test_two_1[[idx]][, 250:1000])))),
  "No pQ (12.01.25)" = colMeans(t(sapply(1:1, function(idx) rowMeans(test_two_2[[idx]][, 250:1000])))),
  .id = "Model"
)

ggplot(
  bind_rows(
    "With pQ (06.01.25)" = do.call(
      rbind, 
      lapply(1:10, function(idx) 
        data.frame(t(test_pQ_1[[idx]])) %>%
          mutate(Time = 1:length(A))
      )
    ),
    "No pQ (06.01.25)" = do.call(
      rbind, 
      lapply(1:10, function(idx) 
        data.frame(t(test_two_1[[idx]])) %>%
          mutate(Time = 1:length(A))
      )
    ),
    "With pQ (12.01.25)" = do.call(
      rbind, 
      lapply(1:10, function(idx) 
        data.frame(t(test_pQ_2[[idx]])) %>%
          mutate(Time = 1:length(A))
      )
    ),
    "No pQ (12.01.25)" = do.call(
      rbind, 
      lapply(1:10, function(idx) 
        data.frame(t(test_two_2[[idx]])) %>%
          mutate(Time = 1:length(A))
      )
    ),
    .id = "Model"
  )  %>%
    filter(Time > 250) %>%
    pivot_longer(cols = c(!Time & !Model)) %>%
    group_by(Time, Model, name) %>% 
    summarise(value = mean(value), .groups = "drop"),
  aes(Time, value, group = Model, color = Model)) +
  geom_line() +
  facet_wrap(vars(name), scales = "free_y")
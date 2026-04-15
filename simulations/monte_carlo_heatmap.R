library(parallel)

# --- Detect cores ---
n_cores <- max(1, detectCores() - 1)
cat("Using", n_cores, "cores\n")

# --- Helper: single simulation run (self-contained for parallel export) ---
run_one_sim_parallel <- function(N, K_true, omega_true,
                                 target_recall = 0.95,
                                 gamma = 0.99,
                                 seed = NULL) {
  
  # Simulate sequence
  if (!is.null(seed)) set.seed(seed)
  
  R_rem <- K_true
  I_rem <- N - K_true
  y <- integer(N)
  
  for (t in seq_len(N)) {
    if (R_rem == 0) { y[t:N] <- 0L; break }
    if (I_rem == 0) { y[t:N] <- 1L; break }
    p_rel <- (omega_true * R_rem) / (omega_true * R_rem + I_rem)
    y[t] <- rbinom(1L, 1L, p_rel)
    if (y[t] == 1L) R_rem <- R_rem - 1L else I_rem <- I_rem - 1L
  }
  
  k_cum <- cumsum(y)
  
  # Coarse grids for speed
  K_grid     <- seq(1, N, by = 20)
  omega_grid <- seq(1, 6, length.out = 15)
  
  # Pre-build parameter matrices (vectorized approach from your other script)
  K_mat     <- matrix(K_grid,     nrow = length(K_grid),
                      ncol = length(omega_grid))
  omega_mat <- matrix(omega_grid, nrow = length(K_grid),
                      ncol = length(omega_grid), byrow = TRUE)
  
  # Uniform joint prior
  post <- matrix(1.0, nrow = length(K_grid), ncol = length(omega_grid))
  post <- post / sum(post)
  
  stop_n <- NA
  
  for (n in seq_len(N)) {
    
    k   <- k_cum[n]
    y_n <- y[n]
    
    # Remaining counts under all (K, omega) — vectorized
    R_rem_mat <- K_mat - (k - y_n)
    I_rem_mat <- (N - K_mat) - ((n - 1L) - (k - y_n))
    
    R_rem_mat[R_rem_mat < 0] <- 0
    I_rem_mat[I_rem_mat < 0] <- 0
    
    denom <- omega_mat * R_rem_mat + I_rem_mat
    p_rel_mat <- ifelse(denom > 0,
                        omega_mat * R_rem_mat / denom,
                        0)
    
    # Likelihood
    lik_mat <- if (y_n == 1L) p_rel_mat else 1 - p_rel_mat
    
    # Update posterior
    post <- post * lik_mat
    s    <- sum(post)
    if (s == 0 || is.na(s)) next
    post <- post / s
    
    # Marginalise over omega to get posterior over K
    post_K <- rowSums(post)
    
    # Posterior probability of meeting recall target
    recall_vals <- k / K_grid
    recall_vals[recall_vals > 1] <- 1
    prob_target <- sum(post_K[recall_vals >= target_recall])
    
    # Stopping rule
    if (!is.na(prob_target) && prob_target >= gamma) {
      stop_n <- n
      break
    }
  }
  
  list(
    recall_at_stop = if (!is.na(stop_n)) k_cum[stop_n] / K_true else NA_real_,
    work_saved     = if (!is.na(stop_n)) 1 - stop_n / N       else 0
  )
}

# --- Experimental setups ---
setups <- list(
  list(name = "Baseline",   N = 10000, prev = 0.05, omega = 3),
  list(name = "Vary omega", N = 10000, prev = 0.05, omega = 5),
  list(name = "Vary prev",  N = 10000, prev = 0.01, omega = 3),
  list(name = "Vary N",     N = 5000,  prev = 0.05, omega = 3)
)

recall_targets <- c(0.90, 0.95, 0.99)
gammas         <- c(0.90, 0.95, 0.99)
n_MC           <- 100

# --- Set up cluster ---
cl <- makeCluster(n_cores)
clusterExport(cl, "run_one_sim_parallel")
clusterEvalQ(cl, set.seed(NULL))

# --- Run experiments ---
results <- list()

for (s in setups) {
  
  N     <- s$N
  K     <- round(s$prev * N)
  omega <- s$omega
  
  cat("Running:", s$name, "\n")
  
  for (tau in recall_targets) {
    for (g in gammas) {
      
      key_name <- paste(s$name, "tau=", tau, "gamma=", g)
      
      # Export current parameters to cluster
      clusterExport(cl,
                    c("N", "K", "omega", "tau", "g"),
                    envir = environment())
      
      # Run all MC iterations in parallel
      mc_results <- parLapply(
        cl,
        seq_len(n_MC),
        function(m) {
          run_one_sim_parallel(
            N           = N,
            K_true      = K,
            omega_true  = omega,
            target_recall = tau,
            gamma       = g,
            seed        = m
          )
        }
      )
      
      results[[key_name]] <- list(
        work_saved = sapply(mc_results, function(x) x$work_saved),
        recall     = sapply(mc_results, function(x) x$recall_at_stop)
      )
      
      cat("  Done:", key_name, "\n")
    }
  }
}

stopCluster(cl)
cat("\nAll experiments complete.\n")

#Adding percentile to table 
for (name in names(results)) {
  
  ws <- results[[name]]$work_saved
  rc <- results[[name]]$recall
  
  cat("\n---", name, "---\n")
  cat("Stopping frequency:", mean(!is.na(rc)), "\n")
  cat("Median work saved:", median(ws, na.rm = TRUE), "\n")
  cat("Median recall:", median(rc, na.rm = TRUE), "\n")
  cat("1% percentile recall:", quantile(rc, 0.01, na.rm = TRUE), "\n")
}


#Building summary dataframe from results list
summary_df <- data.frame()

for (key_name in names(results)) {
  
  ws <- results[[key_name]]$work_saved
  rc <- results[[key_name]]$recall
  
  # Parse the key name - format is "Baseline tau= 0.9 gamma= 0.9"
  # Extract setup name, tau and gamma
  tau_val   <- as.numeric(sub(".*tau= ([0-9.]+).*",   "\\1", key_name))
  gamma_val <- as.numeric(sub(".*gamma= ([0-9.]+).*", "\\1", key_name))
  
  # Extract setup name (everything before " tau=")
  setup_name <- trimws(sub(" tau=.*", "", key_name))
  
  # Map setup name to readable config label
  config_label <- switch(setup_name,
                         "Baseline"   = "N=10000, prev=5%, ω=3 (Baseline)",
                         "Vary omega" = "N=10000, prev=5%, ω=5",
                         "Vary prev"  = "N=10000, prev=1%, ω=3",
                         "Vary N"     = "N=5000,  prev=5%, ω=3",
                         setup_name   # fallback
  )
  
  summary_df <- rbind(summary_df, data.frame(
    config_label      = config_label,
    tau               = tau_val,
    gamma             = gamma_val,
    median_recall     = median(rc, na.rm = TRUE),
    recall_1pct       = quantile(rc, 0.01, na.rm = TRUE),
    median_work_saved = median(ws, na.rm = TRUE),
    miss_rate         = mean(rc < tau_val | is.na(rc), na.rm = FALSE),
    stringsAsFactors  = FALSE
  ))
}

# --- Set factor levels for consistent ordering ---
config_order <- c(
  "N=10000, prev=5%, ω=3 (Baseline)",
  "N=10000, prev=5%, ω=5",
  "N=10000, prev=1%, ω=3",
  "N=5000,  prev=5%, ω=3"
)
summary_df$config_label <- factor(summary_df$config_label,
                                  levels = config_order)

summary_df$tau_gamma_label <- factor(
  paste0("τ=", summary_df$tau, "\nγ=", summary_df$gamma),
  levels = c(
    "τ=0.9\nγ=0.9",  "τ=0.9\nγ=0.95",  "τ=0.9\nγ=0.99",
    "τ=0.95\nγ=0.9", "τ=0.95\nγ=0.95", "τ=0.95\nγ=0.99",
    "τ=0.99\nγ=0.9", "τ=0.99\nγ=0.95", "τ=0.99\nγ=0.99"
  )
)

# --- Plot 1: Heatmap of median recall ---
p1 <- ggplot(summary_df,
             aes(x = tau_gamma_label,
                 y = config_label,
                 fill = median_recall)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", median_recall)),
            size = 3.2, colour = "black") +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = 0.95,
    limits   = c(0.88, 1.00),
    name     = "Median\nRecall"
  ) +
  labs(
    title = "Median Recall at Stopping across Configurations",
    x     = "τ / γ Combination",
    y     = "Experimental Configuration"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x    = element_text(size = 9),
    axis.text.y    = element_text(size = 9),
    plot.title     = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave("plot_mc_recall_heatmap.png", p1,
       width = 10, height = 4, dpi = 150)
cat("Saved: plot_mc_recall_heatmap.png\n")

# --- Plot 2: Heatmap of median work saved ---
p2 <- ggplot(summary_df,
             aes(x = tau_gamma_label,
                 y = config_label,
                 fill = median_work_saved)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", median_work_saved)),
            size = 3.2, colour = "black") +
  scale_fill_gradient2(
    low      = "#d73027",
    mid      = "#ffffbf",
    high     = "#1a9850",
    midpoint = 0.20,
    limits   = c(0.00, 0.55),
    name     = "Median\nWork Saved"
  ) +
  labs(
    title = "Median Work Saved at Stopping across Configurations",
    x     = "τ / γ Combination",
    y     = "Experimental Configuration"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x    = element_text(size = 9),
    axis.text.y    = element_text(size = 9),
    plot.title     = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

ggsave("plot_mc_worksaved_heatmap.png", p2,
       width = 10, height = 4, dpi = 150)
cat("Saved: plot_mc_worksaved_heatmap.png\n")

# --- Plot 3: Heatmap of miss rate ---
p3 <- ggplot(summary_df,
             aes(x = tau_gamma_label,
                 y = config_label,
                 fill = miss_rate)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.3f", miss_rate)),
            size = 3.2, colour = "black") +
  scale_fill_gradientn(
    colours = c("#1a9850", "#ffffbf", "#d73027", "#808080"),
    values  = scales::rescale(c(0, 0.05, 0.15, 1.0)),
    limits  = c(0, 1),
    name    = "Miss\nRate"
  ) +
  labs(
    title = "Miss Rate across Configurations",
    x     = "τ / γ Combination",
    y     = "Experimental Configuration"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x    = element_text(size = 9),
    axis.text.y    = element_text(size = 9),
    plot.title     = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# Flag cells where stopping never triggered
summary_df <- summary_df %>%
  mutate(label_text = ifelse(miss_rate == 1.000,
                             "N/A\n(never\nstopped)",
                             sprintf("%.3f", miss_rate)))

# Then in the geom_text line, replace the label argument:
# aes(label = sprintf("%.3f", miss_rate))
# with:
# aes(label = label_text)

ggsave("plot_mc_missrate_heatmap.png", p3,
       width = 10, height = 4, dpi = 150)
cat("Saved: plot_mc_missrate_heatmap.png\n")

# --- Print summary table to console ---
cat("\nFull summary table:\n")
print(
  summary_df %>%
    select(config_label, tau, gamma,
           median_recall, recall_1pct,
           median_work_saved, miss_rate) %>%
    arrange(config_label, tau, gamma),
  row.names = FALSE
)


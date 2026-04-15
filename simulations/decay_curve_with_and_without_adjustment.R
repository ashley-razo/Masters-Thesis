#function defining the bias behavior
omega_true_function <- function(n, N, remaining_rel, K_true, scenario, decay_strength = 1) {
  
  if (scenario == "constant") {
    omega <- 3
  }
  
  if (scenario == "linear_decay") {
    omega_start <- 5
    omega_end <- 5 - decay_strength * 3 #controls strength of decay
    omega <- omega_start - (omega_start - omega_end) * (n/N)
  }
  
  if (scenario == "sudden_drop") {
    if (n < 0.3*N) {
      omega <- 5
    } else {
      omega <- 2
    }
  }
  
  if (scenario == "remaining_rel") {
    omega <- 1 + 4 * (remaining_rel / K_true)
  }
  
  return(omega)
}

#biased urn generator
simulate_sequence <- function(N, K_true, scenario, 
                              decay_strength = 1, seed=NULL){
  
  if(!is.null(seed)) set.seed(seed)
  
  remaining_rel <- K_true
  remaining_irrel <- N - K_true
  
  y <- numeric(N)
  
  for(n in 1:N){
    
    if(remaining_rel == 0){
      y[n] <- 0
      remaining_irrel <- remaining_irrel - 1
      next
    }
    
    if(remaining_irrel == 0){
      y[n] <- 1
      remaining_rel <- remaining_rel - 1
      next
    }
    
    omega_current <- omega_true_function(
      n, N, remaining_rel, K_true, scenario, decay_strength
    )
    
    p_rel <- (omega_current * remaining_rel) /
      (omega_current * remaining_rel + remaining_irrel)
    
    y[n] <- rbinom(1,1,p_rel)
    
    if(y[n]==1){
      remaining_rel <- remaining_rel - 1
    } else{
      remaining_irrel <- remaining_irrel - 1
    }
    
  }
  
  return(y)
}


#Bayesian inference and stopping rule 
run_one_sim <- function(N, K_true,
                        target_recall = 0.95,
                        gamma = 0.99,
                        scenario,
                        decay_strength = 1,
                        bias_factor = 1,
                        seed = NULL) {
  
  y <- simulate_sequence(N, K_true, scenario, decay_strength, seed)
  k_cum <- cumsum(y)
  
  
  K_grid <- seq(1, N, by = 5)
  
  omega_grid <- seq(1, 6, length.out = 40)
  
  # Joint prior (uniform)
  post <- matrix(1, nrow = length(K_grid), ncol = length(omega_grid))
  post <- post / sum(post)
  
  stop_n <- NA
  
  for (n in 1:N) {
    
    k_prev <- if (n == 1) 0 else k_cum[n-1]
    n_prev <- n - 1
    
    R_rem <- outer(K_grid - k_prev,
                   rep(1, length(omega_grid)))
    R_rem[R_rem < 0] <- 0
    
    I_rem <- outer(N - K_grid - (n_prev - k_prev),
                   rep(1, length(omega_grid)))
    I_rem[I_rem < 0] <- 0
    
    omega_mat <- matrix(rep(omega_grid,
                            each = length(K_grid)),
                        nrow = length(K_grid))

    omega_mat <- bias_factor * omega_mat
    
    
    
    p_rel <- (omega_mat * R_rem) /
      (omega_mat * R_rem + I_rem)
    
    p_rel[is.na(p_rel)] <- 0
    
    if (y[n] == 1) {
      post <- post * p_rel
    } else {
      post <- post * (1 - p_rel)
    }
    
    post <- post / sum(post)
    
    # posterior over K
    post_K <- rowSums(post)
    
    recall_probs <- k_cum[n] / K_grid
    recall_probs[recall_probs > 1] <- 1
    
    prob_target <- sum(post_K[recall_probs >= target_recall])
    
    if (prob_target >= gamma) {
      stop_n <- n
      break
    }
  }
  
  recall_at_stop <- if (!is.na(stop_n))
    k_cum[stop_n] / K_true else NA
  
  work_saved <- if (!is.na(stop_n))
    1 - stop_n / N else 0
  
  return(list(
    recall_at_stop = recall_at_stop,
    work_saved = work_saved
  ))
}

#experimental setups
scenarios <- c(
  "constant",
  "linear_decay",
  "sudden_drop",
  "remaining_rel"
)

#monte carlo for each setup
n_MC <- 100

N <- 10000
prev <- 0.05
K <- round(prev * N)

decay_strengths <- c(0.25, 0.5, 0.75, 1)

results_decay <- list()

for (ds in decay_strengths) {
  
  cat("Running decay strength:", ds, "\n")
  
  recall_vec <- numeric(n_MC)
  
  for (m in 1:n_MC) {
    
    out <- run_one_sim(
      N = N,
      K_true = K,
      scenario = "linear_decay",
      decay_strength = ds,
      bias_factor = 0.8,
      seed = m
    )
    
    recall_vec[m] <- out$recall_at_stop
  }
  
  results_decay[[as.character(ds)]] <- recall_vec
}

#extracting 1% recall 
recall_1pct_values <- sapply(results_decay, function(x) quantile(x, 0.01))



# Parameters
N <- 10000
prev <- 0.05
K <- round(prev * N)
n_MC <- 100
decay_strengths <- c(0.25, 0.50, 0.75, 1.00)
target_recall <- 0.95
gamma <- 0.99

#Running WITHOUT conservative adjustment (bias_factor = 1.0)
recall_1pct_unadjusted <- numeric(length(decay_strengths))

for (i in seq_along(decay_strengths)) {
  ds <- decay_strengths[i]
  cat("Unadjusted | Decay strength:", ds, "\n")
  recall_vec <- numeric(n_MC)
  
  for (m in 1:n_MC) {
    out <- run_one_sim(
      N = N,
      K_true = K,
      scenario = "linear_decay",
      decay_strength = ds,
      bias_factor = 1.0,
      seed = m,
      target_recall = target_recall,
      gamma = gamma
    )
    recall_vec[m] <- out$recall_at_stop
  }
  recall_1pct_unadjusted[i] <- quantile(recall_vec, 0.01)
}

#Running WITH conservative adjustment (bias_factor = 0.8)
recall_1pct_adjusted <- numeric(length(decay_strengths))

for (i in seq_along(decay_strengths)) {
  ds <- decay_strengths[i]
  cat("Adjusted | Decay strength:", ds, "\n")
  recall_vec <- numeric(n_MC)
  
  for (m in 1:n_MC) {
    out <- run_one_sim(
      N = N,
      K_true = K,
      scenario = "linear_decay",
      decay_strength = ds,
      bias_factor = 0.8,
      seed = m,
      target_recall = target_recall,
      gamma = gamma
    )
    recall_vec[m] <- out$recall_at_stop
  }
  recall_1pct_adjusted[i] <- quantile(recall_vec, 0.01)
}

# --- Plot both curves on the same axes ---
png("plot_decay_comparison.png", width = 700, height = 500, res = 100)

plot(
  decay_strengths, recall_1pct_unadjusted,
  type = "b",
  pch = 16,
  col = "black",
  lwd = 2,
  ylim = c(0.90, 0.97),
  xlab = "Decay Strength",
  ylab = "1% Recall at Stopping",
  main = "Effect of Conservative Adjustment on Recall Guarantee\nunder Linear Bias Decay"
)

lines(
  decay_strengths, recall_1pct_adjusted,
  type = "b",
  pch = 17,
  col = "steelblue",
  lwd = 2,
  lty = 2
)

abline(h = 0.95, col = "red", lwd = 2, lty = 1)

legend(
  "bottomleft",
  legend = c(
    "Without adjustment (bias factor = 1.0)",
    "With adjustment (bias factor = 0.8)",
    "Recall target (τ = 0.95)"
  ),
  col = c("black", "steelblue", "red"),
  lwd = c(2, 2, 2),
  lty = c(1, 2, 1),
  pch = c(16, 17, NA),
  bty = "n"
)

dev.off()

cat("\nSummary of results:\n")
cat(sprintf("%-12s %-20s %-20s\n",
            "Decay", "1% Recall (Unadj)", "1% Recall (Adj)"))
for (i in seq_along(decay_strengths)) {
  cat(sprintf("%-12.2f %-20.4f %-20.4f\n",
              decay_strengths[i],
              recall_1pct_unadjusted[i],
              recall_1pct_adjusted[i]))
}
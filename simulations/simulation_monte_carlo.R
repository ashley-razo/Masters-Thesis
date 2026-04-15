#biased urn generator
simulate_sequence <- function(N, K_true, omega_true, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  remaining_rel  <- K_true
  remaining_irrel <- N - K_true
  
  y <- numeric(N)
  
  for (n in 1:N) {
    
    if (remaining_rel == 0) {
      y[n] <- 0
      remaining_irrel <- remaining_irrel - 1
      next
    }
    
    if (remaining_irrel == 0) {
      y[n] <- 1
      remaining_rel <- remaining_rel - 1
      next
    }
    
    p_rel <- (omega_true * remaining_rel) /
      (omega_true * remaining_rel + remaining_irrel)
    
    y[n] <- rbinom(1, 1, p_rel)
    
    if (y[n] == 1) {
      remaining_rel <- remaining_rel - 1
    } else {
      remaining_irrel <- remaining_irrel - 1
    }
  }
  
  return(y)
}

#Bayesian inference and stopping rule 
run_one_sim <- function(N, K_true, omega_true,
                        target_recall = 0.95,
                        gamma = 0.99,
                        seed = NULL,
                        plot = FALSE) {
  
  y <- simulate_sequence(N, K_true, omega_true, seed)
  k_cum <- cumsum(y)
  
  # Parameter grids
  #K_grid <- 1:N        #taking too long 
  
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
    stop_n = stop_n,
    recall_at_stop = recall_at_stop,
    work_saved = work_saved
  ))
}

#experimental setups
setups <- list(
  list(name = "Baseline", N = 10000, prev = 0.05, omega = 3),
  list(name = "Vary omega", N = 10000, prev = 0.05, omega = 5),
  list(name = "Vary prev", N = 10000, prev = 0.01, omega = 3),
  list(name = "Vary N", N = 5000, prev = 0.05, omega = 3)
)

#monte carlo for each setup
n_MC <- 100   # increase if runtime allows

#NEW: vary stopping parameters
recall_targets <- c(0.9, 0.95, 0.99)
gammas <- c(0.9, 0.95, 0.99)
###
results <- list()

for (s in setups) {
  
  cat("Running:", s$name, "\n")
  
  N <- s$N
  K <- round(s$prev * N)
  omega <- s$omega
  
  #NEW: for varying parameters
  for (tau in recall_targets) {
    for (g in gammas) {
      key_name <- paste(s$name,
                        "tau=", tau,
                        "gamma=", g)
      
      work_saved_vec <- numeric(n_MC)
      recall_vec <- numeric(n_MC)
      
      for (m in 1:n_MC) {
        
        out <- run_one_sim(
          N, K, omega,
          target_recall = tau,
          gamma = g,
          seed = m
        )
        
        work_saved_vec[m] <- out$work_saved
        recall_vec[m] <- out$recall_at_stop
      }
      
      results[[key_name]] <- list(
        work_saved = work_saved_vec,
        recall = recall_vec
      )
      
      cat("  Done:", key_name, "\n")
    }
  }
}
  
  
  #work_saved_vec <- numeric(n_MC)
  #recall_vec <- numeric(n_MC)
  
  #for (m in 1:n_MC) {
    
   # out <- run_one_sim(N, K, omega, seed = m)
    
    #work_saved_vec[m] <- out$work_saved
    #recall_vec[m] <- out$recall_at_stop
  #}
  
  #results[[s$name]] <- list(
   # work_saved = work_saved_vec,
    #recall = recall_vec
  #)
#}


#plot distributions
for (name in names(results)) {
  
  ws <- results[[name]]$work_saved
  rc <- results[[name]]$recall
  
  # Work saved distribution
  hist(ws,
       main = paste("Work Saved -", name),
       xlab = "Work Saved",
       breaks = 20)
  
  # Recall distribution
  hist(rc,
       main = paste("Recall at Stop -", name),
       xlab = "Recall",
       breaks = 20,
       xlim = c(0.9, 1))
  abline(v = 0.95, col = "red", lwd = 2)
}

#Adding percentile to table 
for (name in names(results)) {
  
  ws <- results[[name]]$work_saved
  rc <- results[[name]]$recall
  
  cat("\n---", name, "---\n")
  cat("Stopping frequency:", mean(rc > 0), "\n")
  cat("Median work saved:", median(ws), "\n")
  cat("Median recall:", median(rc), "\n")
  cat("1% percentile recall:", quantile(rc, 0.01), "\n")
}



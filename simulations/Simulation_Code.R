rm(list = ls())

#simulating a biased urn draw sequence (without replacement)
simulate_biased_urn <- function(N, K, omega = 1, n_draws = N, seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  stopifnot(N >= 1, K >= 0, K <= N, omega > 0)
  n_draws <- min(n_draws, N)
  
  R_rem <- K
  I_rem <- N - K
  
  y <- integer(n_draws) # 1 = relevant, 0 = irrelevant 
  
  for (t in seq_len(n_draws)) {
    #If one class exhausted, remaining draws are forced 
    if (R_rem == 0) {
      y[t:n_draws] <- 0
      break
    }
    if (I_rem == 0) {
      y[t:n_draws] <- 1
      break
    }
    
    p_rel <- (omega * R_rem) / (omega * R_rem + I_rem)
    y[t] <- rbinom(1, size = 1, prob = p_rel)
    
    if (y[t] == 1) R_rem <- R_rem - 1 else I_rem <- I_rem - 1
  }
  
  y
} 

#Bayesian posterior for K in the unbiased case (w = 1)
posterior_K_unbiased <- function(N, n_screened, k_found, prior = NULL) {
  #Support constraints:
  #K must be at least k_found 
  #K cannot exceed k_found + (N - n_screened) because only N - n_screened remain unseen 
  K_min <- max(1, k_found) #make sure k is at least 1 (so never dividing by zero when computing recall)
  K_max <- k_found + (N - n_screened)
  K_grid <- K_min:K_max
  
  #using a uniform distribution as prior for K
  if (is.null(prior)) {
    prior <- rep(1 / length(K_grid), length(K_grid)) #uniform on allowed support 
  } else {
    stopifnot(length(prior) == length(K_grid))
    prior <- prior / sum(prior)
  }
  
  lik <- dhyper(x = k_found, m = K_grid, n = N - K_grid, k = n_screened)
  post_unnorm <- lik * prior
  post <- post_unnorm / sum(post_unnorm)
  
  list(K_grid = K_grid, posterior = post)
}

#Credible interval for recall and posterior "safe to stop" probability 
credible_interval_recall <- function(K_grid, posterior, k_found, level = 0.95) {
  recall_vals <- k_found / K_grid
  #recall decreases as K increases; reorder by recall for quantiles 
  o <- order(recall_vals)
  recall_vals <- recall_vals[o]
  posterior <- posterior[o]
  
  cdf <- cumsum(posterior)
  alpha <- (1 - level) / 2
  
  lo <- recall_vals[which(cdf >= alpha) [1]]
  hi <- recall_vals[which(cdf >= 1 - alpha) [1]]
  
  c(lower = lo, upper = hi)
}

prob_recall_ge <- function(K_grid, posterior, k_found, target_recall) {
  recall_vals <- k_found / K_grid
  sum(posterior[recall_vals >= target_recall])
}


#Parameters 
N <- 2000
K_true <- 120
omega_true <- 1   #unbiased case for the "unit test"
target_recall <- 0.95
stop_conf <- 0.99 #require 99% posterior probability of meeting target 

#Simulate Screening 
y <- simulate_biased_urn(N, K_true, omega = omega_true, n_draws = N, seed = 1)
  #implementing monte carlo simulation of biased urn model
k_cum <- cumsum(y) #k_cum[n] shows number of relevant documents found after screening n records 
n_seq <- seq_along(y) #n_seq is the screening step index(1,2, ... , N)

#Run Bayesian updating over time 
stop_n <- NA
ci_store <- matrix(NA, nrow = N, ncol = 2)
p_safe_store <- rep(NA, N)
K_map_store <- rep(NA, N)

for (n in n_seq) {
  k <- k_cum[n]
  
  #using Bayesian updating to estimate posterior distribution of K
  post <- posterior_K_unbiased(N, n_screened = n, k_found = k)
  K_grid <- post$K_grid
  pK <- post$posterior
  
  #MAP estimate for K (most likely K)
  K_map <- K_grid[which.max(pK)]
  K_map_store[n] <- K_map
  
  #Credible interval for recall
  ci <- credible_interval_recall(K_grid, pK, k_found = k, level = 0.95)
  ci_store[n, ] <- ci
  
  #Probability that recall >= target 
  p_safe <- prob_recall_ge(K_grid, pK, k_found = k, target_recall = target_recall)
  p_safe_store[n] <- p_safe
  
  if(is.na(stop_n) && !is.na(p_safe) && p_safe >= stop_conf) {
    stop_n <- n
  }
}

cat("True K =", K_true, "\n")
cat("Stop n =", stop_n, "\n")
cat("At stop: k_found =", k_cum[stop_n], 
    "empirical recall =", k_cum[stop_n]/K_true, "\n")
cat("Posterior P(recall >= target) =", p_safe_store[stop_n], "\n")

#Simple plots
#Plot that checks as screening progresses, posterior should peak around true K
plot(n_seq, K_map_store, type = "l",
     xlab = "Documents screened (n)", ylab = "MAP estimate of K",
     main = "Posterior MAP of total relevant K over screening")
abline(h = K_true, lty = 2)
if (!is.na(stop_n)) abline(v = stop_n, lty = 3)

plot(n_seq, p_safe_store, type = "l",
     xlab = "Documents screened (n)", 
     ylab = paste0("P(recall ≥ ", target_recall, ")"),
     main = "Posterior safety probability over screening")
abline(h = stop_conf, lty = 2)
if (!is.na(stop_n)) abline(v = stop_n, lty = 3)

plot(n_seq, ci_store[,1], type = "l",
     xlab = "Documents screened (n)", ylab = "Recall",
     main = "95% credible interval for recall")
lines(n_seq, ci_store[,2], type = "l")
abline(h = target_recall, lty = 2)

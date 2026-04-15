ls()
#rm(list = ls())
library("BiasedUrn")

#simulating a biased urn draw sequence (without replacement)
#helper functions:
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
  if (all(is.na(recall_vals))) return(NA_real_)
  sum(posterior[recall_vals >= target_recall], na.rm = TRUE)
}

#Define credible interval for K
credible_interval_K <- function(K_grid, posterior, level = 0.95) {
  cdf <- cumsum(posterior)
  alpha <- (1 - level) / 2
  
  lo <- K_grid[which(cdf >= alpha)[1]]
  hi <- K_grid[which(cdf >= 1 - alpha)[1]]
  
  c(lower = lo, upper = hi)
}

#Define credible interval for ω 
credible_interval_omega <- function(omega_grid, posterior, level = 0.95) {
  cdf <- cumsum(posterior) 
  alpha <- (1 - level) / 2
  
  lo <- omega_grid[which(cdf >= alpha) [1]]
  hi <- omega_grid[which(cdf >= 1 - alpha) [1]]
  
  c(lower = lo, upper = hi)
}

#Simulation engine
run_one_sim <- function(N, K_true, omega_true) {

#Parameters 
target_recall <- 0.95
stop_conf <- 0.99 #require 99% posterior probability of meeting target 

#Simulate Screening 
y <- simulate_biased_urn(N, K_true, omega = omega_true, seed = 1)
#implementing monte carlo simulation of biased urn model
k_cum <- cumsum(y) #k_cum[n] shows number of relevant documents found after screening n records 
n_seq <- seq_along(y) #n_seq is the screening step index(1,2, ... , N)

#Run Bayesian updating over time 
#Storage vectors
stop_n <- NA
p_safe_store <- rep(NA, N)
K_map_store <- rep(NA, N)
omega_map_store <- rep(NA, N)
ciRecall95 <- ciRecall80 <- matrix(NA, N, 2)
ciK95 <- ciK80  <- matrix(NA, N, 2)
ciOmega95 <- ciOmega80 <- matrix(NA, N, 2)

# Fixed grids for entire run
K_grid <- 0:N
omega_grid <- seq(1, 5, length.out = 50)

# Prior (only once)
prior_K     <- rep(1/length(K_grid), length(K_grid))
prior_omega <- rep(1/length(omega_grid), length(omega_grid))
post_joint  <- outer(prior_K, prior_omega)


#Bayesian Loop
#for (n in n_seq[n_seq %% 10 == 0]) {
for (n in n_seq) {
  k <- k_cum[n]
  
  #compute likelihood
  lik_mat <- matrix(NA, nrow = length(K_grid), ncol = length(omega_grid))
  
  y_n <- y[n]   # current draw (0 or 1)
  
  lik_mat <- matrix(0, length(K_grid), length(omega_grid))
  
  for (i in seq_along(K_grid)) {
    for (j in seq_along(omega_grid)) {
      
      K  <- K_grid[i]
      ω  <- omega_grid[j]
      
      if (K < k || K > N - (n - k)) next
      
      R_rem <- K - k + y_n
      I_rem <- (N - K) - ((n-1) - (k - y_n))
      
      if (R_rem < 0 || I_rem < 0) next
      
      p_rel <- (ω * R_rem) / (ω * R_rem + I_rem)
      
      lik_mat[i,j] <- ifelse(y_n == 1, p_rel, 1 - p_rel)
    }
  }
  
  # Sequential Bayes update
  post_joint <- post_joint * lik_mat
  
  if (sum(post_joint) == 0 || is.na(sum(post_joint))) next
  post_joint <- post_joint / sum(post_joint)
  
  
  #marginalize
  post_K <- rowSums(post_joint)
  ciK95[n, ] <- credible_interval_K(K_grid, post_K, 0.95)
  ciK80[n, ] <- credible_interval_K(K_grid, post_K, 0.80)
  post_omega <- colSums(post_joint)
  
  #MAP estimate of omega
  omega_map_store[n] <- omega_grid[which.max(post_omega)]
  
  #Credible interval for omega
  ciOmega95[n, ] <- credible_interval_omega(omega_grid, post_omega, 0.95)
  ciOmega80[n, ] <- credible_interval_omega(omega_grid, post_omega, 0.80)
  
  # MAP estimate of K
  K_map_store[n] <- K_grid[which.max(post_K)]
  
  # Recall credible interval
  ciRecall95[n, ] <- credible_interval_recall(K_grid, post_K, k, 0.95)
  ciRecall80[n, ] <- credible_interval_recall(K_grid, post_K, k, 0.80)
  
  # Posterior safety probability
  p_safe_store[n] <- prob_recall_ge(
    K_grid, post_K, k_found = k, target_recall
  )
  
  #Stopping rule
  if (is.na(stop_n) &&
      !is.na(p_safe_store[n]) &&
      p_safe_store[n] >= stop_conf) {
    stop_n <- n
  }
}
  

cat("True K =", K_true, "\n")
cat("Stop n =", stop_n, "\n")
cat("At stop: k_found =", k_cum[stop_n], 
    "empirical recall =", k_cum[stop_n]/K_true, "\n")
cat("Posterior P(recall >= target) =", p_safe_store[stop_n], "\n")

#Simple plots

par(mfrow = c(3,1), mar = c(4,4,2,1))

# K plot
plot(n_seq, K_map_store, type="l",
     ylim = range(ciK95, ciK80, K_true, na.rm=TRUE),
     ylab="Total relevant K", xlab="Documents screened",
     main="Posterior K")

# 95% CI
lines(n_seq, ciK95[,1], lty=2)
lines(n_seq, ciK95[,2], lty=2)

# 80% CI
lines(n_seq, ciK80[,1], lty=3)
lines(n_seq, ciK80[,2], lty=3)

# True K
abline(h = K_true, col=2, lty=3)

# stopping time
if(!is.na(stop_n)) abline(v=stop_n, lty=3)


#Recall plot
recall_map <- k_cum / K_map_store
recall_true <- k_cum / K_true

plot(n_seq, recall_map, type="l",
     ylim = c(0,1),
     ylab="Recall", xlab="Documents screened",
     main="Recall")

# 95% CI
lines(n_seq, ciRecall95[,1], lty=2)
lines(n_seq, ciRecall95[,2], lty=2)

# 80% CI
lines(n_seq, ciRecall80[,1], lty=3)
lines(n_seq, ciRecall80[,2], lty=3)

# True recall
lines(n_seq, recall_true, col=2)

# target
abline(h = target_recall, lty=3)

# stopping
if(!is.na(stop_n)) abline(v=stop_n, lty=3)


#Omega plot
plot(n_seq, omega_map_store, type="l",
     ylim = range(ciOmega95, ciOmega80, omega_true, na.rm=TRUE),
     ylab="Bias ω", xlab="Documents screened",
     main="Bias parameter")

lines(n_seq, ciOmega95[,1], lty=2)
lines(n_seq, ciOmega95[,2], lty=2)

lines(n_seq, ciOmega80[,1], lty=3)
lines(n_seq, ciOmega80[,2], lty=3)

abline(h = omega_true, col=2, lty=3)
if(!is.na(stop_n)) abline(v=stop_n, lty=3)

} 

#end of simulation engine 

#Defining realistic dataset scenarios 
scenarios <- list(
  list(N = 1000, K = 50, omega = 2),
  list(N = 2000, K = 120, omega = 3),
  list(N = 5000, K = 250, omega = 4),
  list(N = 10000, K = 500, omega = 5)
)

#running all scenarios 
par(mfrow = c(length(scenarios), 3), mar = c(4,4,2,1))

for(s in scenarios) {
  run_one_sim(s$N, s$K, s$omega)
}

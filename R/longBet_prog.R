#' Demo model with time-varying treatment effect with prognostic effect
#' 
#' Fit panel data with individual random effect, time random effect and 
#' heterogeneous treatment effect with XBCF.
#' y_it = alpha_i + gamma_t + mu(x, t) + tau(x_i, k)*z_it + epsilon_it
#' 
#' @param y Numeric matrix of observed outcomes
#' @param x Numeric matrix of static covariates
#' @param z Numeric vector of treatment variables
#' @param t0 Numeric value of treatment time
#' 
#' @return An object of fitted parameters
#' 
#' @examples 
#' 
#' @export
longBet_prog <- function(y, x, z, t0, mc = 100, burnin = 10, ntrees = 10,
                       mu_a = 0, nu_a = 1, alpha_a = 3, beta_a = 10, 
                       mu_g = 0, nu_g = 1, alpha_g = 3, beta_g = 10,
                       a = 16, b = 4){
  
  # # prior params
  # mu_a = 0; nu_a = 1; alpha_a = 3; beta_a = 10
  # mu_g = 0; nu_g = 1; alpha_g = 3; beta_g = 10
  # a = 16; b = 4
  
  # required package
  # library(XBCF)
  
  t1 <- ncol(y)
  
  # # standardize per unit pre-treatment period
  # mean_y = rowMeans(y[,1:(t0-1)])
  # sd_y = apply(y[,1:(t0-1)], 1, sd)
  # for (i in 1:n){
  #   y[i,] = (y[i,] - mean_y[i])/sd_y[i]
  # }
  
  # ini parameters
  tauhat <- array(0, dim = c(mc, n, t1))
  muhat <- array(0, dim = c(mc, n, t1))
  sigma2hat <- rep(0, mc)
  
  tau_vec <- rep(0,n*t1)
  tau_mat <- matrix(tau_vec, n, t1)
  mu_vec <- rep(0, n*t1)
  mu_mat <- matrix(mu_vec, n, t1)
  
  xbcf_y <- as.vector(y)
  mean_y <- mean(xbcf_y)
  x_con <- c()
  x_mod <- c()
  for (i in 1:t1){
    x_con <- rbind(x_con, cbind(x, i))
    x_mod <- rbind(x_mod, cbind(x, max(i - (t0 - 1), 0)))
  }
  xbcf_z <- c(rep(0, n*(t0-1)), rep(z, t1-t0+1))
  xbcf.fit <- XBCF(as.matrix(xbcf_y), as.matrix(xbcf_z), as.matrix(x_con), as.matrix(x_mod),
                   num_sweeps = mc, burnin = burnin, n_trees_con = ntrees, n_trees_mod = ntrees,
                   pcat_con = 0, pcat_mod = 0)
  
  tauhats.adjusted <- xbcf.fit$tauhats[,(burnin+1):mc] * xbcf.fit$b_draws[(burnin+1):mc, 2]
  muhats.adjusted <- xbcf.fit$muhats.adjusted - xbcf.fit$tauhats[,(burnin+1):mc] * xbcf.fit$b_draws[(burnin+1):mc, 1]
    
  # tauhats.adjusted <- xbcf.fit$tauhats.adjusted
  # muhats.adjusted <- xbcf.fit$muhats.adjusted 
  
  for (i in 1:(mc-burnin)){
    tauhat[i,,] <- matrix(tauhats.adjusted[,i], n, t1)
    muhat[i,,] <- matrix(xbcf.fit$muhats.adjusted[,i], n, t1)
  }
    
  obj <- list()
  obj$tauhat <- tauhat
  obj$muhat <- muhat
  
  return(obj)
  
}
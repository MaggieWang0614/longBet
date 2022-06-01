#' Demo model with time-varying treatment effect using post-treatment period as covariate
#' 
#' Fit panel data with XBCF
#' y_it = mu(x, t) + tau(x_i, k)*z_it + epsilon_it
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
longBet_xbcf <- function(y, x, z, t0, mc = 100, burnin = 10, ntrees = 10, 
                         alpha_mod = 0.95, beta_mod = 1.25){

  # required package
  # library(XBCF)
  
  t1 <- ncol(y)
  tauhat <- array(0, dim = c(mc - burnin, n, t1))
  muhat <- array(0, dim = c(mc - burnin, n, t1))
  
  xbcf_y <- as.vector(y)
  mean_y <- mean(xbcf_y)
  x_con <- c()
  x_mod <- c()
  for (i in 1:t1){
    x_con <- rbind(x_con, cbind(x, i))
    x_mod <- rbind(x_mod, cbind(x, max(0, i - t0 + 1)))
  }
  xbcf_z <- c(rep(0, (t0-1)*n), rep(z, t1-t0+1))
  xbcf.fit <- XBCF(as.matrix(xbcf_y), as.matrix(xbcf_z), as.matrix(x_con), as.matrix(x_mod),
                   num_sweeps = mc, burnin = burnin, n_trees_con = 0, n_trees_mod = ntrees,
                   pcat_con = 0, pcat_mod = 0, alpha_mod = alpha_mod, beta_mod = beta_mod)
  tauhats.adjusted <- xbcf.fit$tauhats[, (burnin+1):mc] * xbcf.fit$b_draws[(burnin+1):mc, 2]
  for (i in 1:(mc-burnin)){
    tauhat[i,,] <- matrix(tauhats.adjusted[,i], n, t1) + mean_y
    muhat[i,,] <- matrix(xbcf.fit$muhats.adjusted[,i], n, t1)
  }
  
  obj <- list()
  obj$tauhat <- tauhat
  obj$muhat <- muhat
  
  return(obj)
  
}
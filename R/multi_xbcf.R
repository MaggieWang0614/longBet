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
multi_xbcf <- function(y, x, z, t0, mc = 100, burnin = 10, ntrees = 10){
  
  # required package
  # library(XBCF)
  
  t1 <- ncol(y)
  tauhat <- array(0, dim = c(mc - burnin, n, t1 - t0 + 1))
  muhat <- array(0, dim = c(mc - burnin, n, t1 - t0 + 1))
  
  for (t in t0:t1){
    xbcf_y <- as.matrix(y[,t])
    mean_y <- mean(xbcf_y)
    xbcf_x <- as.matrix(x)
    xbcf_z <- as.matrix(z) 
    
    xbcf.fit <- XBCF(as.matrix(xbcf_y), as.matrix(xbcf_z), as.matrix(xbcf_x), as.matrix(xbcf_x),
                     num_sweeps = mc, burnin = burnin, n_trees_con = 0, n_trees_mod = ntrees,
                     pcat_con = 0, pcat_mod = 0, alpha_mod = 0.95, beta_mod = 1.25)
  
    tauhats.adjusted <- xbcf.fit$tauhats[, (burnin + 1):mc] * xbcf.fit$b_draws[(burnin + 1):mc,2]
    tauhat[,,t - t0 + 1] <- t(tauhats.adjusted)  + mean_y
    muhat[,,t - t0 + 1] <- t(xbcf.fit$muhats.adjusted)
  }
  
  obj <- list()
  obj$tauhat <- tauhat
  obj$muhat <- muhat
  
  return(obj)
  
}
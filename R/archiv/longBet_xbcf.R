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
longBet_xbcf <- function(y, x, z, t0, mc = 100, burnin = 10, ntrees = 10){

  # required package
  # library(XBCF)
  
  t1 <- ncol(y)
  
  xbcf_y <- as.vector(y)
  x_con <- c()
  x_mod <- c()
  for (i in 1:t1){
    x_con <- rbind(x_con, cbind(x, i))
    x_mod <- rbind(x_mod, cbind(x, max(0, i - t0 + 1)))
  }
  xbcf_z <- c(rep(0, (t0-1)*n), rep(z, t1-t0+1))
  xbcf.fit <- XBCF(as.matrix(xbcf_y), as.matrix(xbcf_z), as.matrix(x_con), as.matrix(x_mod),
                   num_sweeps = mc, burnin = burnin, n_trees_con = ntrees, n_trees_mod = ntrees,
                   pcat_con = 0, pcat_mod = 0)
  
  tauhat <- array(0, dim = c(n, t1, mc))
  muhat <- array(0, dim = c(n, t1, mc))
  for (i in 1:mc){
    tauhat[,,i] <- matrix(xbcf.fit$tauhats[,i], n, t1)
    muhat[,,i] <- matrix(xbcf.fit$muhats[,i], n, t1)
  }
  
  obj <- list()
  obj$tauhat <- tauhat[, , (burnin+1):mc]
  obj$muhat <- muhat[, ,(burnin+1):mc]
  
  return(obj)
  
}
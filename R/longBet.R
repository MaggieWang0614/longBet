#' Demo model with time-varying treatment effect using post-treatment period as covariate
#' 
#' Fit panel data with individual random effect, time random effect and 
#' heterogeneous treatment effect with XBCF.
#' y_it = alpha_i + gamma_t + tau(x_i, k)*z_it + epsilon_it
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
longBet_fixte <- function(y, x, z, t0, mc = 100, burnin = 10, ntrees = 10,
                          mu_a = 0, nu_a = 1, alpha_a = 3, beta_a = 2, 
                          mu_g = 0, nu_g = 1, alpha_g = 3, beta_g = 2,
                          a = 16, b = 4){
  
  # prior params
  mu_a = 0; nu_a = 1; alpha_a = 3; beta_a = 10
  mu_g = 0; nu_g = 1; alpha_g = 3; beta_g = 10
  a = 16; b = 4
  
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
  alphahat <- matrix(0, n, mc)
  gammahat <- matrix(0, t1, mc)
  tauhat <- array(0, dim = c(mc, n, t1-t0+1))
  sigma2hat <- rep(0, mc)
  
  alpha_vec <- rep(0, n)
  gamma_vec <- rep(0, t1)
  tau_vec <- rep(0,n*(t1-t0+1))
  tau_mat <- matrix(tau_vec, n, t1-t0+1)
  
  # hierachical params
  mu_alpha <- sigma_alpha2 <- mu_gamma <- sigma_gamma2 <- rep(0, mc)
  
  # ini params
  sigma2 <- 1 / rgamma(1, a, b)
  
  res <- y # as if res = y - 0 - 0 - 0 * z
  
  # iteration
  for (iter in 1:mc){
    # update hierachical param for alpha
    # alpha_i ~ N(mu_alpha, sigma_alpha)
    # (mu_alpha, sigma_alpha2) ~ Normal-IG(mu_a, nu_a, alpha_a, beta_a)
    # sigma_alpha2  | alpha_i ~ posterior IG
    post_alpha_a <- alpha_a + n / 2
    post_beta_a <- beta_a + 0.5*sum((alpha_vec - mean(alpha_vec))^2) + 0.5*n*nu_a * (mean(alpha_vec) - mu_a)^2 / (n + nu_a)
    sigma_alpha2[iter] <- 1 / rgamma(1, post_alpha_a, post_beta_a)
    # mu_alpha | sigma_alpha2, mu_a, nu_a ~ N(mu_a, sigma_alpha2/nu_a)
    post_mu_a <- (nu_a * mu_a + n*mean(alpha_vec)) / (nu_a + n)
    post_nu_a <-  nu_a + n 
    mu_alpha[iter] <- rnorm(1, post_mu_a, sqrt(sigma_alpha2[iter]/post_nu_a))
    
    # update hierachical param for gamma
    # gamma_t ~ N(mu_gamma, sigma_gamma)
    # (mu_gamma, sigma_gamma2) ~ Normal-IG(mu_g, nu_g, alpha_g, beta_g)
    # sigma_gamma2  | gamma_t ~ posterior IG
    post_alpha_g <- alpha_g + t1 / 2
    post_beta_g <- beta_g + 0.5*sum((gamma_vec - mean(gamma_vec))^2) + 0.5*t1*nu_g * (mean(gamma_vec) - mu_g)^2 / (t1 + nu_g)
    sigma_gamma2[iter] <- 1 / rgamma(1, post_alpha_g, post_beta_g)
    # mu_gamma | sigma_gamma2, mu_g, nu_g ~ N(mu_g, sigma_gamma2/nu_g)
    post_mu_g <- (nu_g * mu_g + t1*mean(gamma_vec)) / (nu_g + t1)
    post_nu_g <-  nu_g + t1
    mu_gamma[iter] <- rnorm(1, post_mu_g, sqrt(sigma_gamma2[iter]/post_nu_g))
    
    # update alpha
    res <- res + matrix(rep(alpha_vec, t1), n, t1) # update residual, add alpha from prev iter
    for (i in 1:n){
      # alpha_i ~ N(mu_alpha, sigma_alpha2)
      # res | alpha_i ~ N(alpha_i, sigma2)
      # alpha_i | res ~ N() # as if conjugate normal with known variance
      post_mu_alpha <- (mu_alpha[iter] / sigma_alpha2[iter] + sum(res[i,]) / sigma2) / (1/sigma_alpha2[iter] + n/sigma2)
      post_sigma_alpha2 <- 1 / (1 / sigma_alpha2[iter] + n / sigma2)
      alpha_vec[i] <- rnorm(1, post_mu_alpha, sqrt(post_sigma_alpha2))
    }
    res <- res - matrix(rep(alpha_vec, t1), n, t1) # update residual
    alphahat[, iter] <- alpha_vec
    
    # update gamma
    res <- res + t(matrix(rep(gamma_vec, n), t1, n))
    for (i in 1:t1){
      # gamma_t ~ N(mu_gamma, sigma_gamma2)
      # res | gamma_t ~ N(gamma_t, sigma2)
      post_mu_gamma <- (mu_gamma[iter] / sigma_gamma2[iter] + sum(res[,i]) / sigma2) / (1/sigma_gamma2[iter] + n/sigma2)
      post_sigma_gamma2 <- 1 / (1 / sigma_gamma2[iter] + n / sigma2)
      gamma_vec[i] <- rnorm(1, post_mu_gamma, sqrt(post_sigma_gamma2))
    }
    res <- res - t(matrix(rep(gamma_vec, n), t1, n)) # update residual
    gammahat[, iter] <- gamma_vec
    
    # update treatment effect
    # a demo model, input tau(x, k) for y_{i, t0+k}
    # vectorize residuals
    res[,t0:t1] <- res[, t0:t1] + tau_mat*matrix(rep(z, t1-t0+1), n, t1-t0+1)
    xbcf_y <- as.vector(res[,t0:t1])
    xbcf_x <- c()
    for (i in 1:(t1-t0+1)){
      xbcf_x <- rbind(xbcf_x, cbind(x, i))
    }
    xbcf_z <- rep(z, t1-t0+1)
    xbcf.fit <- XBCF(as.matrix(xbcf_y), as.matrix(xbcf_z), as.matrix(xbcf_x), as.matrix(xbcf_x),
                    num_sweeps = 1, burnin = 0, n_trees_con = 0, n_trees_mod = ntrees,
                    pcat_con = 0, pcat_mod = 0)
    
    tau_vec <- xbcf.fit$tauhats
    tau_mat <- matrix(tau_vec, n, t1-t0+1)
    res[,t0:t1] <- res[, t0:t1] - tau_mat*matrix(rep(z, t1-t0+1), n, t1-t0+1)
    tauhat[iter,,] <- tau_mat
    
    
    # update sigma
    # sigma2 ~ IG(a, b)
    # res | sigma2 ~ N(0, sigma2)
    # sigma2 | res ~ IG(a + n/2, b + sum((res - 0)^2) / 2)
    sigma2 <- 1 / rgamma(1, a + n*t1/2, b + sum(res^2)/2)
    sigma2hat[iter] <- sigma2
  }
  
  obj = list()
  obj$alphahat <- alphahat[, (burnin+1):mc]
  obj$gammahat <- gammahat[, (burnin+1):mc]
  obj$tauhat <- tauhat[(burnin+1):mc, , ]
  obj$sigma2hat <- sigma2hat[(burnin+1):mc]
  
  return(obj)
  
}
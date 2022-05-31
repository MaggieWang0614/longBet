setwd('~/Dropbox (ASU)/longBet/R/')
## ----------------------------------------------------------------------------------------------------------------
n = 100
t1 = 40
t0 = 21
p = 1
mu_a = 0
sig_a = 1
mu_g = 1
sig_g = 1


mc <- 100
tau_true_mc <- array(0, dim = c(mc, n, t1-t0+1))
tau_hat_mc <- array(0, dim = c(mc, n, t1-t0+1))
pct_bias_mc <- array(0, dim = c(mc, n, t1-t0+1))

for (iter in 1:100){
  alpha = rnorm(n, mu_a, sig_a)
  gamma = rnorm(t1, mu_g, sig_g)
  ## ----------------------------------------------------------------------------------------------------------------
  x = as.matrix(rnorm(n))
  tau <- function(x){a = (x + 1.5)^2; 5*sqrt(a) + sin(5*a) + 1}
  tau_mat <- matrix(0, n, t1-t0+1)
  tau_mat[,1] <- tau(x)
  for (i in 2:(t1-t0+1)){
    tau_mat[,i] <- 0.9*tau_mat[,i-1]
  }
  ## ----------------------------------------------------------------------------------------------------------------
  z = rbinom(n, 1, 0.5)
  
  
  ## ----------------------------------------------------------------------------------------------------------------
  eps = matrix(rnorm(n*t1, 0, 0.2), nrow = n, ncol = t1)
  
  
  ## ----------------------------------------------------------------------------------------------------------------
  y0 = y1 = y = matrix(0, nrow = n, ncol = t1)
  for (i in 1:n){
    y0[i,] = y0[i,] + alpha[i]
  }
  for (j in 1:t1){
    y0[,j] = y0[,j] + gamma[j]
  }
  y0 = y0 + eps
  y1 = y0
  y1[, t0:t1] = y0[, t0:t1] + tau_mat
  z_mat = matrix(rep(z, t1), n, t1)
  y = y0 * (1-z_mat) + y1 * z_mat
  
  
  ## ----------------------------------------------------------------------------------------------------------------
  source('longBet_re.R')
  library(XBCF)
  fit <- longBet_re(y, x, z, t0, 50, 20)
  
  
  ## ----------------------------------------------------------------------------------------------------------------
  tau_hat <- colMeans(fit$tauhat)
  pct_bias <- abs((tau_hat - tau_mat) / tau_mat)
  tau_true_mc[iter,,] <- tau_mat
  tau_hat_mc[iter,,] <- tau_hat
  pct_bias_mc[iter,,] <- pct_bias
  
  save(tau_true_mc, tau_hat_mc, pct_bias_mc, file = "longBet_sim.RData")
  
}

png("Percentage bias on 100 draws.png", width = 600, height = 300)
pct_bias <- colMeans(pct_bias_mc)
par(mfrow=c(1,2))
plot(t0:t1, colMeans(pct_bias), type = "l", col = 1, ylim = range(0, colMeans(pct_bias), 1),lty = 2, ylab = 'Percentage bias on tau over time', yaxt="n")
# lines(t0:t1, rep(0, t1-t0+1), col = 3, lty = 1, lwd = 2)
axis(2, at=pretty(colMeans(pct_bias)), lab=pretty(colMeans(pct_bias)) * 100, las=TRUE)
legend("topleft", legend = c("pct bias"), col = c(1), lty = c(1,2))

plot(x, rowMeans(pct_bias), col = 1, ylim = range(0, rowMeans(pct_bias), 1), ylab = 'Percentage bias on tau over individual', yaxt="n")
axis(2, at=pretty(rowMeans(pct_bias)), lab=pretty(rowMeans(pct_bias)) * 100, las=TRUE)
# legend("topleft", legend = c("pct_bias"), col = c(1), lty = c(1,2))
dev.off()


# Generate dgp (residuals)
t1 <- 100 # time
n <- 1000 # observation per time
sig <- 0.5
lambda <- 2
sig2 <- sig^2
lambda2 <- lambda^2

# generate time coefficient
time_coef <- function(t) {100 * dpois(t + 1, lambda = 20)}
sqrd_exp <- function(x1, x2) {
  sig2*exp( -(x1 - x2)^2 / (2*lambda2) )
}

t <- 1:t1
mu_t <- time_coef(t) # mean
Sigma_t <- outer(t, t, sqrd_exp)

L <- svd(Sigma_t)$u %*% diag(sqrt(svd(Sigma_t)$d))
beta <- mu_t + L %*% rnorm(t)

# generate noise and observations
sig_param <- rgamma(n*t, shape = 0.1*sd(beta), rate = 1) # variance of noise
eps <- rnorm(n*t, 0, sig_param) # noise
y <- rep(beta, n) + eps

# GP with basis expansion
time_b <- proc.time()
k <- t1/2 # number of knots / basis
mu <- seq(1, t1, length.out = k) # centers of basis
sigma_w <- 1 # prior variance for w ~ N(0, sigma_w^2*I)
Lambda0 <- solve(sigma_w^2 * diag(rep(1,k)))

sqrd_exp_b <- function(t){
  (pi / 2)^(-1/4) * lambda^(-1/2) * exp(- (t - mu)^2 / lambda2)
}
Sigma_mu <- t(sapply(t, sqrd_exp_b))

# Bayseian linear regression
X <- do.call(rbind, replicate(n, Sigma_mu, simplify=FALSE))
Lambda_N <- t(X) %*% X + Lambda0
inv_LambdaN <- solve(Lambda_N)
mu_N <- inv_LambdaN %*% (t(X) %*% y)
L_N <- svd(inv_LambdaN)$u %*% diag(sqrt(svd(inv_LambdaN)$d))
w <- mu_N + L_N %*% rnorm(k)
beta_blr <- Sigma_mu %*% w
time_b <- proc.time() - time_b
cat('rmse basis:', sqrt(mean((beta_blr - beta)^2)), '\n')
cat('time basis:', time_b[1])

## Benchmark with GP
time_gp <- proc.time()
sigma <- mean(sig_param)
inv_Sigma <- solve(Sigma_t)
Sigma_gp <- solve(inv_Sigma + n * sigma^(-2) * diag(rep(1, nrow(inv_Sigma))))
mu_gp <- Sigma_gp %*% (n * sigma^(-2) * rowMeans(matrix(y, t1, n)))
svd_gp <- svd(Sigma_gp)
L_gp <- svd_gp$u %*% diag(sqrt(svd_gp$d))
beta_gp <- mu_gp + L_gp %*% rnorm(t1)
time_gp <- proc.time() - time_gp
cat('rmse gp:', sqrt(mean((beta_gp - beta)^2)),'\n')
cat('time gp:', time_gp[1])

rm(list=ls())
v<-data.frame(df_mfgarch[,1:11])

## function for MIDAS GARCH estimation 

# Phi
calculate_phi <- function(w1, w2, K) {
  weights <- sapply(c(1:K),
                    FUN = function(j) (j / (K + 1))^(w1 - 1) * (1 - j / (K + 1))^(w2 - 1))
  weights <- weights/sum(weights)
  weights
}

# estimation of tau_t
# c++ function
sum_tau <- function(m, theta, phivar, covariate, K) {
  .Call('_mfGARCH_sum_tau', PACKAGE = 'mfGARCH', m, theta, phivar, covariate, K)
}

# estimation of tau
calculate_tau_mf <- function(df, x, low.freq, w1, w2, theta, m, K){
  phi.var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), x)
  tau <- c(rep(NA, times = K),
           exp(sum_tau(m = m, theta = theta, phivar = phi.var, covariate = x, K = K)))
  
  result <- merge(df, cbind(unique(df[low.freq]), tau), by = low.freq)
  result
}

# estimation g_i,t
# c++ function
calculate_g <- function(omega, alpha, beta, gamma, returns, g0) {
  .Call('_mfGARCH_calculate_g', PACKAGE = 'mfGARCH', omega, alpha, beta, gamma, returns, g0)
}


## components estimation using log-liklelyhood ##
# Normal
llh_simple_g <- function(y, mu, alpha, beta, gamma, m, g_zero) {
  # parametri per funzione
  omega <- 1 - alpha - beta - gamma / 2
  ret <- y
  ret_std <- (ret - mu)/sqrt(exp(m))
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ret_std, g0 = g_zero)
  # funzione di log-verosimiglianza della normale
  ll<-1/2 * log(2 * pi) + 1/2 * log(g * exp(m)) + 1/2 * (ret - mu)^2/(g * exp(m))
  # return della funzione
  return(ll)
}

# 
lfg <- function(p) {
  # get the log-likelihood
  llf <- llh_simple_g(y = v$return,
                      mu     = p["mu"],
                      alpha  = p["alpha"],
                      beta   = p["beta"],
                      gamma  = p["gamma"],
                      m      = p["m"],
                      g_zero = var(v$return))
  # return output
  return(llf)
}

par.start <- c(mu = 0.00, alpha = 0.01, beta = 0.85, gamma = 0.04, m = 0)
ui.opt    <- rbind(c(0, -1, -1, -1/2, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
ci.opt    <- c(-0.999, 0, 0)
maxit     <- 100

p.e.nlminb_g <- constrOptim(theta = par.start,
                          f = function(theta) { sum(lfg(theta)) },
                          control = list(maxit = maxit, reltol = 1e-4, trace = TRUE),
                          grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)


# t-Student
llh_simple_t <- function(y, mu, alpha, beta, gamma, m, g_zero) {
  # gradi di liberta
  dgrs<-30.0
  # componenti 
  omega <- 1 - alpha - beta - gamma / 2
  ret <- y
  ret_std <- (ret - mu)/sqrt(exp(m))
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ret_std, g0 = g_zero)
  # funzione di verosimiglianza
  #llt<- 
       - gammaln_internal(0.5 * (dgrs + 1.0)) + gammaln_internal(0.5 * dgrs) 
       + 0.5 * log(g * exp(m) * pi * (dgrs - 2.0)) 
       + (0.5 * (dgrs + 1.0)) * log(1.0 + ((ret - mu)^2 / (g * exp(m) * (dgrs - 2.0))))
  # return(llt)
  # non parte se uso return
}

lft <- function(p) {
  
  # get the log-likelihood
  llf <- llh_simple_t(y = v$return,
                      mu     = p["mu"],
                      alpha  = p["alpha"],
                      beta   = p["beta"],
                      gamma  = p["gamma"],
                      m      = p["m"],
                      g_zero = var(v$return))
  return(llf)
}

par.start <- c(mu = 0.00, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0.00)
ui.opt    <- rbind(c(0, -1, -1, -1/2, 0), c(0, 1, 0, 0, 0), 
                   c(0, 0, 1, 0, 0), c(0, 0, 0, 1, 0))
ci.opt    <- c(-0.999, 0, 0, 0)
maxit     <- 100

p.e.nlminb <- constrOptim(theta = par.start,
                          f = function(theta) { sum(lft(theta)) },
                          control = list(maxit = maxit, reltol = 1e-4, trace = TRUE),
                          grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)






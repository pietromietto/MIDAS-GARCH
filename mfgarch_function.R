rm(list=ls())
## MIDAS GARCH model estimation
## Assuming Gaussian likelyhood

# components estimation

# phi estimation
calculate_phi <- function(w1, w2, K) {
  weights <- sapply(c(1:K),
                    FUN = function(j) (j / (K + 1))^(w1 - 1) * (1 - j / (K + 1))^(w2 - 1))
  weights <- weights/sum(weights)
  weights
}

# tau_t estimation
sum_tau <- function(m, theta, phivar, covariate, K) {
  .Call('_mfGARCH_sum_tau', PACKAGE = 'mfGARCH', m, theta, phivar, covariate, K)
}

calculate_tau_mf <- function(df, x, low.freq, w1, w2, theta, m, K){
  phi.var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), x)
  tau <- c(rep(NA, times = K),
           exp(sum_tau(m = m, theta = theta, phivar = phi.var, covariate = x, K = K)))
  
  result <- merge(df, cbind(unique(df[low.freq]), tau), by = low.freq)
  result
}

sum_tau(m = 0.1, theta = 0.6, phivar = 0.4, covariate = v$dhousing, K = 1)


# g_i,t estimation
calculate_g <- function(omega, alpha, beta, gamma, returns, g0) {
  .Call('_mfGARCH_calculate_g', PACKAGE = 'mfGARCH', omega, alpha, beta, gamma, returns, g0)
}


# simple case
llh_simple <- function(y, mu, alpha, beta, gamma, m, g_zero) {
  omega <- 1 - alpha - beta - gamma / 2
  ret <- y
  ret_std <- (ret - mu)/sqrt(exp(m))
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ret_std, g0 = g_zero)
  1/2 * log(2 * pi) + 1/2 * log(g * exp(m)) + 1/2 * (ret - mu)^2/(g * exp(m))
}

llh_simple<-function(y, mu, alpha, beta, gamma, m, g_zero){
  omega<-1-alpha-beta-gamma/2
  ret<-y
  ret_std<-(ret-mu)/sqrt(exp(m))
  g<-calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                 returns = ret_std, g0 = g_zero)
  ll<-gammaln((4+1)/2)-gammaln(4/2)- 0.5*log(g*exp(m)*pi*(4-2)) - ((4+1)/2) * log(1+((ret-mu)^2/(g*exp(m)*(4-2))))
  return(ll)
}

lf <- function(p) {
  llh_simple(y = v$return,
             mu = p["mu"],
             alpha = p["alpha"],
             beta = p["beta"],
             gamma = p["gamma"],
             m = p["m"],
             g_zero = var(v$return))
}
par.start <- c(mu = 0.00, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0)
ui.opt <- rbind(c(0.00, -1, -1, -1/2, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
ci.opt <- c(-0.99999, 0, 0)

p.e.nlminb <- constrOptim(theta = par.start,
                          f = function(theta) { sum(lf(theta)) },
                          grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)

par <- -p.e.nlminb$par
round(par, 3)

# case with tau_t
llh_mf <-function(df, x, y, low.freq, mu, omega, alpha, beta, gamma,
           m, theta, w1 = 1, w2 = 1, g_zero, K){
             
    tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
    
    ret <- y
    ret <- ret[which.min(is.na(tau)):length(ret)]  # lags can't be used for likelihood
    tau <- tau[which.min(is.na(tau)):length(tau)]
    
    g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                     returns = ((ret - mu)/sqrt(tau)), g0 = g_zero)
    
    if (sum(g <= 0) > 0) {
      #rep(NA, times = length(y))
      #stop("g_t seems to be negative for at least one point in time?")
      rep(NA, times = length(g))
    } else {
      1/2 * log(2 * pi) + 1/2 * log(g * tau) + 1/2 * (ret - mu)^2/(g * tau)
    }
  }

# K==1
df_llh <- v[, c("return", "dhousing", "year_month")]
df_llh[, "year_month"] <- as.integer(unlist(df_llh[ , "year_month"]))

lf <- function(p) {
  llh_mf(df = df_llh, y = "return", x = "dhousing", low.freq = "year_month",
         mu = p["mu"],
         omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
         alpha = p["alpha"],
         beta = p["beta"],
         gamma = p["gamma"],
         m = p["m"],
         theta = p["theta"],
         w1 = 1, w2 = 1, g_zero = var(v$return), K = 1)
}
par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0)
ui_opt <- rbind(c(0, -1, -1, -1/2, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
ci_opt <- c(-0.99999, 0, 0)

p.e.nlminb <- constrOptim(theta = par_start, f = function(theta) { sum(lf(theta)) },
                          grad = NULL, ui = ui_opt, ci = ci_opt, hessian = FALSE)  

# 
low_freq_backup <- data[, low.freq]
if (x != "date") {
  if (is.null(x.two) == TRUE) {
    df_llh <- data[, c(y, x, low.freq)]
    df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
  } else {
    low_freq.two_backup <- data[, low.freq.two]
    if (low.freq != low.freq.two) { # if they are different, both have to be included in df_llh
      df_llh <- data[, c(y, x, low.freq, x.two, low.freq.two)]
      df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
      df_llh[, low.freq.two] <- as.integer(unlist(df_llh[ , low.freq.two]))
    } else { # else, the low.freq column is needed only once
      df_llh <- data[, c(y, x, low.freq, x.two)]
      df_llh[, low.freq] <- as.integer(unlist(df_llh[ , low.freq]))
    }
  }
}

df_llh <- v[, c("return", "dhousing", "year_month")]
df_llh[, "year_month"] <- as.integer(unlist(df_llh[ , "year_month"]))
df_llh[, "dhousing"]<-as.double(df_llh[, "dhousing"])

str(df_llh$dhousing)

gamma(5)

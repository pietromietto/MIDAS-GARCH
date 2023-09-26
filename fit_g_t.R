rm(list=ls())
# 
library(numDeriv)

### FILE R ###
### Upload c++ files
calculate_g <- function(omega, alpha, beta, gamma, returns, g0) {
  .Call('_mfGARCH_calculate_g', PACKAGE = 'mfGARCH', omega, alpha, beta, gamma, returns, g0)
}

calculate_h_andersen <- function(ndays, delta, mu, theta, omega, lambda, Z, pi, h0) {
  .Call('_mfGARCH_calculate_h_andersen', PACKAGE = 'mfGARCH', ndays, delta, mu, theta, omega, lambda, Z, pi, h0)
}

calculate_p <- function(ndays, delta, mu, Zp, h, p0) {
  .Call('_mfGARCH_calculate_p', PACKAGE = 'mfGARCH', ndays, delta, mu, Zp, h, p0)
}

simulate_r <- function(n_days, n_intraday, alpha, beta, gamma, Z, h0) {
  .Call('_mfGARCH_simulate_r', PACKAGE = 'mfGARCH', n_days, n_intraday, alpha, beta, gamma, Z, h0)
}

simulate_r_rv_as_dependent <- function(n_days, n_intraday, alpha, beta, gamma, Z, h0, K, m, theta, weights, lowfreq, rvol) {
  .Call('_mfGARCH_simulate_r_rv_as_dependent', PACKAGE = 'mfGARCH', n_days, n_intraday, alpha, beta, gamma, Z, h0, K, m, theta, weights, lowfreq, rvol)
}

sum_tau_fcts <- function(i, m, theta, phivar, covariate, K) {
  .Call('_mfGARCH_sum_tau_fcts', PACKAGE = 'mfGARCH', i, m, theta, phivar, covariate, K)
}

sum_tau <- function(m, theta, phivar, covariate, K) {
  .Call('_mfGARCH_sum_tau', PACKAGE = 'mfGARCH', m, theta, phivar, covariate, K)
}

### Carico funzioni utili per il calcolo delle componenti
forecast_garch <- function(omega, alpha, beta, gamma, g, ret, steps.ahead) {
  omega / (1 - alpha - gamma/2 - beta) + (alpha + beta + gamma/2)^(steps.ahead - 1) * (omega + (alpha + gamma/2 * as.numeric(ret < 0)) * ret^2 + beta * g - omega / (1 - alpha - gamma/2 - beta))
}

calculate_phi <- function(w1, w2, K) {
  weights <- sapply(c(1:K),
                    FUN = function(j) (j / (K + 1))^(w1 - 1) * (1 - j / (K + 1))^(w2 - 1))
  weights <- weights/sum(weights)
  weights
}

calculate_tau <- function(covariate, w1, w2, theta, m, K) { # used for simulation
  phi_var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), covariate)
  tau <- c(rep(NA, times = K),
           exp(sum_tau(m = m, theta = theta, phivar = phi_var, covariate = covariate, K = K)))
  tau
}

calculate_tau_mf <- function(df, x, low.freq, w1, w2, theta, m, K,
                             x.two = NULL, K.two = NULL, theta.two = NULL,
                             low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
  phi.var <- calculate_phi(w1, w2, K)
  covariate <- c(rep(NA, times = K), x)
  tau <- c(rep(NA, times = K),
           exp(sum_tau(m = m, theta = theta, phivar = phi.var, covariate = x, K = K)))
  
  result <- merge(df, cbind(unique(df[low.freq]), tau), by = low.freq)
  
  if (is.null(x.two) == FALSE) {
    phi.var.two <- calculate_phi(w1.two, w2.two, K.two)
    covariate.two <- c(rep(NA, times = K.two), x.two)
    tau.two <- c(rep(NA, times = K.two),
                 exp(sum_tau(m = 0, theta = theta.two, phivar = phi.var.two,
                             covariate = x.two, K = K.two)))
    result <- merge(result, cbind(unique(df[low.freq.two]), tau.two), by = low.freq.two)
    
    result$tau.one <- result$tau # store tau component due to first covariate
    result$tau <- result$tau.one * result$tau.two # generate joint tau component
  }
  
  result
}

### CASO NORMALE ###
### (funzioni che poi sono da modificare con inserimento della distribuziobne t-Student)
llh_mf_g <-
  function(df, x, y, low.freq, mu, omega, alpha, beta, gamma,
           m, theta, w1 = 1, w2 = 1, g_zero, K,
           x.two = NULL, K.two = NULL, theta.two = NULL,
           low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
    
    if (is.null(x.two) == FALSE) {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K,
                              x.two = x.two, K.two = K.two, theta.two = theta.two,
                              low.freq.two = low.freq.two,
                              w1.two = w1.two, w2.two = w2.two)$tau
    } else {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
    }
    
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

llh_simple_g <- function(y, mu, alpha, beta, gamma, m, g_zero) {
  omega <- 1 - alpha - beta - gamma / 2
  ret <- y
  ret_std <- (ret - mu)/sqrt(exp(m))
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ret_std, g0 = g_zero)
  1/2 * log(2 * pi) + 1/2 * log(g * exp(m)) + 1/2 * (ret - mu)^2/(g * exp(m))
}

### FUNZIONE fit_mfGARCH

fit_mfgarch_g <- function(data, y, x = NULL, K = NULL, low.freq = "date", var.ratio.freq = NULL, gamma = TRUE, weighting = "beta.restricted", x.two = NULL, K.two = NULL, low.freq.two = NULL, weighting.two = NULL, multi.start = FALSE, control = list(par.start = NULL)) {
  
  print("For ensuring numerical stability of the parameter optimization and inversion of the Hessian, it is best to multiply log returns by 100.")
  
  if (is.null(weighting.two) == FALSE) {
    if (weighting.two != "beta.restricted") {
      stop("Right now, only beta.restricted weighting scheme is employed for the second covariate.")
    }
  }
  
  if (is.null(x.two) == FALSE) {
    weighting.two <- "beta.restricted"
  }
  
  if (is.null(x.two) == FALSE && gamma == FALSE) {
    stop("Regarding two covariates, only asymmetric GJR-GARCH component is implemented.")
  }
  
  # if (K == 1 && is.null(x.two) == FALSE && K.two != 1) {
  #   stop("Regarding two covariates, only K.two = 1 is implemented if K = 1.")
  # }
  if (is.null(K.two) == FALSE) {
    if (K == 1 & K.two > 1) {
      stop("Regarding two covariates with one of them being equal to one, only K.two = 1 is implemented.")
    }
  }
  
  if (is.null(x.two) == FALSE) {
    print("Specifying two covariates may lead to long estimation times.")
  }
  
  if (weighting %in% c("beta.restricted", "beta.unrestricted") == FALSE) {
    stop("Incorrect weighting scheme specified - options are \"beta.restricted\" and \"beta.unrestricted\".")
  }
  if (gamma %in% c(TRUE, FALSE) == FALSE) {
    stop("Gamma can't be anything different than TRUE or FALSE.")
  }
  if ("date" %in% colnames(data) == FALSE) {
    stop("No date column.")
  }
  if (inherits(data$date, 'Date') == FALSE) {
    stop(paste0("Supplied date column is not of format 'Date'. It is of class '", class(data$date), "'."))
  }
  if (inherits(data[[low.freq]], 'Date') == FALSE) {
    stop(paste0("Supplied low.freq column is not of format 'Date'. It is of class '", class(data[[low.freq]]), "'."))
  }
  if (is.null(x) == FALSE && K == 0) {
    warning("You specified an external covariate x but chose K = 0 - simple GARCH is estimated (K = 0).")
  }
  
  if (is.null(x) == TRUE) {
    warning("No external covariate x is specified - simple GARCH is estimated (K = 0).")
    x <- "date"
    K <- 0
  }
  if (is.null(K) == TRUE) {
    warning("No K is specified - simple GARCH is estimated (K = 0).")
    x <- "date"
    K <- 0
  }
  if (K < 0 || K %% 1 != 0) {
    stop("K can't be smaller than 0 and has to be an integer.")
  }
  if (dim(unique(data[c(x, low.freq)]))[1] > dim(unique(data[c(low.freq)]))[1]) {
    stop("There is more than one unique observation per low frequency entry.")
  }
  # if ((is.null(x) == TRUE && (is.null(K) == TRUE)) || K == 0) {
  #   K <- 0
  # }
  if (y %in% colnames(data) == FALSE) {
    stop(paste("There is no variable in your data frame with name ", y, "."))
  }
  if (x %in% colnames(data) == FALSE && is.null(x) != FALSE) {
    stop(paste("There is no variable in your data frame with name ", x, "."))
  }
  if (low.freq %in% colnames(data) == FALSE) {
    stop(paste("There is no low freq. variable in your data frame with name ", low.freq, "."))
  }
  if ("tau" %in% colnames(data) == TRUE) {
    stop("There may not be a column named tau - it will be part of df.fitted")
  }
  if ("g" %in% colnames(data) == TRUE) {
    stop("There may not be a column named g - it will be part of df.fitted")
  }
  if (is.null(x) == TRUE) {
    if (sum(is.na(data[[y]]) == TRUE) > 0) {
      stop(paste0("Column ", y, "contains NAs."))
    }
  } else {
    if (sum(is.na(data[[y]]) == TRUE) > 0 | sum(is.na(data[[x]]) == TRUE) > 0) {
      stop(paste0("Either column ", y, " or column ", x, "includes NAs."))
    }
  }
  if (length(unlist(unique(data[["date"]]))) != dim(data)[1]) {
    stop("There is more than one observation per high frequency (presumably date).")
  }
  if (is.null(var.ratio.freq) == FALSE) {
    if (var.ratio.freq %in% colnames(data) == FALSE) {
      stop(paste0("There is no var.ratio.freq column with name ", var.ratio.freq, "."))
    }
  }
  
  # Order by high frequency variable
  data <- data[order(data$date), ]
  # Deprecated dplyr version
  #data <- dplyr::arrange_(data, "date")
  # We store date in new variable because computation on integerized date seemed to be faster
  date_backup <- data[["date"]]
  data["date"] <- as.numeric(unlist(data["date"]))
  
  if (is.null(var.ratio.freq) == TRUE) {
    var.ratio.freq <- low.freq
    print(paste0("No frequency specified for calculating the variance ratio - default: low.freq = ", low.freq))
  }
  
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
  
  g_zero <- var(unlist(data[[y]]))
  ret <- data[[y]]
  
  # Parameter estimation
  if (K == 0) {
    if (gamma == TRUE) {
      lf <- function(p) {
        llh_simple_g(y = ret,
                   mu = p["mu"],
                   alpha = p["alpha"],
                   beta = p["beta"],
                   gamma = p["gamma"],
                   m = p["m"],
                   g_zero = g_zero)
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0)
      ui.opt <- rbind(c(0, -1, -1, -1/2, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
      ci.opt <- c(-0.99999, 0, 0)
    } else {
      lf <- function(p) {
        llh_simple_g(y = ret,
                   mu = p["mu"],
                   alpha = p["alpha"],
                   beta = p["beta"],
                   gamma = 0,
                   m = p["m"],
                   g_zero = g_zero)
        
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0)
      ui.opt <- rbind(c(0, -1, -1, 0), c(0, 1, 0, 0), c(0, 0, 1, 0))
      ci.opt <- c(-0.99999, 0, 0)
    }
    
    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }
    
    p.e.nlminb <- constrOptim(theta = par.start,
                              f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    
    if (multi.start == TRUE && gamma == TRUE) {
      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.two) == "try-error") {
        print("Second-step BFGS optimization failed. Fallback: First-stage Nelder-Mead estimate.")
      } else {
        if (p.e.nlminb.two$value < p.e.nlminb$value) {
          p.e.nlminb <- p.e.nlminb.two
        }
      }
    }
    
    p.e.nlminb$value <- -p.e.nlminb$value
    
    par <- p.e.nlminb$par
    returns <- as.numeric(unlist(data[[y]]))
    tau <- rep(exp(par["m"]), times = length(returns))
    
    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] -  par["gamma"]/2,
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
      tau <- rep(exp(par["m"]), times = length(g))
      
    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = 0,
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
      tau <- rep(exp(par["m"]), times = length(g))
    }
    
    if ((var.ratio.freq %in% c("date", "low.freq")) == FALSE) {
      df.fitted <- cbind(data[c("date", y, var.ratio.freq)], g = g, tau = tau)
    } else {
      df.fitted <- cbind(data[c("date", y)], g = g, tau = tau)
    }
    
    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
  } else { # if K > 0 we get the covariate series
    covariate <- unlist(unique(data[c(low.freq, x)])[x])
    
    if (is.null(x.two) == FALSE) {
      covariate.two <- unlist(unique(data[c(low.freq.two, x.two)])[x.two])
    }
  }
  
  if (K == 1) {
    if (is.null(K.two) == FALSE) {
      if (gamma == TRUE) {
        lf <- function(p) {
          llh_mf_g(df = df_llh, y = ret, x = covariate, low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = p["gamma"],
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K,
                 x.two = covariate.two,
                 K.two = K.two, low.freq.two = low.freq.two,
                 theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0, theta.two = 0)
        ui_opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0), c(0, 1, 0, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      } else {
        lf <- function(p) {
          llh_mf_g(df = df_llh,
                 y = ret, x = covariate, low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = 0,
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K,
                 x.two = covariate.two,
                 K.two = K.two, low.freq.two = low.freq.two,
                 theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, theta.two = 0)
        ui_opt <- rbind(c(0, -1, -1,  0, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      }
      
    } else {
      
      if (gamma == TRUE) {
        lf <- function(p) {
          llh_mf_g(df = df_llh, y = ret, x = covariate, low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = p["gamma"],
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0)
        ui_opt <- rbind(c(0, -1, -1, -1/2, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      } else {
        lf <- function(p) {
          llh_mf_g(df = df_llh,
                 y = ret, x = covariate, low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = 0,
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0)
        ui_opt <- rbind(c(0, -1, -1,  0, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      }
    }
    
    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }
    
    p.e.nlminb <- constrOptim(theta = par_start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui_opt, ci = ci_opt, hessian = FALSE)
    par <- p.e.nlminb$par
    p.e.nlminb$value <- -p.e.nlminb$value
    
    if (is.null(x.two) == FALSE) {
      tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                              w1 = 1, w2 = 1, theta = par["theta"], m = par["m"], K = K,
                              x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                              low.freq.two = low.freq.two,
                              w1.two = 1, w2.two = 1)$tau
    } else {
      tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                              theta = par["theta"], m = par["m"], w1 = 1, w2 = 1, K = K)$tau
    }
    
    tau_forecast <-
      exp(sum_tau_fcts(m = par["m"],
                       i = K + 1,
                       theta = par["theta"],
                       phivar = calculate_phi(w1 = 1, w2 = 1, K = K),
                       covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                       K = K))
    
    if (is.null(x.two) == FALSE) {
      tau_forecast <-
        tau_forecast *
        exp(sum_tau_fcts(m = 0,
                         i = K.two + 1,
                         theta = par["theta.two"],
                         phivar = calculate_phi(w1 = 1, w2 = 1, K = K.two),
                         covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                         K = K.two))
    }
    
    
    returns <- unlist(data[y])
    
    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                         alpha = par["alpha"], beta = par["beta"], gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))), g0 = g_zero))
    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],
                         alpha = par["alpha"], beta = par["beta"], gamma = 0,
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))), g0 = g_zero))
    }
    
    if ((var.ratio.freq %in% c("date", "low.freq")) == FALSE) {
      df.fitted <- cbind(data[c("date", y, low.freq, x, var.ratio.freq)], g = g, tau = tau)
    } else {
      df.fitted <- cbind(data[c("date", y, low.freq, x)], g = g, tau = tau)
    }
    
    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
    
  }
  
  if (K > 1) {
    if (gamma == TRUE) {
      if (weighting == "beta.restricted" & is.null(K.two) == TRUE) {
        lf <- function(p) {
          llh_mf_g(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = p["gamma"],
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1,
                 w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                       m = 0, theta = 0, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0),
                        c(0,  0,  0,    0, 0, 0, 1),
                        c(0,  1,  0,    0, 0, 0, 0),
                        c(0,  0,  1,    0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 0, 0)
      }
      if (weighting == "beta.restricted" & is.null(K.two) == FALSE) {
        if (K.two == 1) {
          lf <- function(p) {
            llh_mf_g(df = df_llh,
                   y = ret,
                   x = covariate,
                   low.freq = low.freq,
                   mu = p["mu"],
                   omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                   alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                   m = p["m"], theta = p["theta"],
                   w1 = 1, w2 = p["w2"], g_zero = g_zero, K = K,
                   x.two = covariate.two,
                   K.two = 1, low.freq.two = low.freq.two,
                   theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
          }
          par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                         m = 0, theta = 0, w2 = 3, theta.two = 0)
          ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 1, 0),
                          c(0,  1,  0,    0, 0, 0, 0, 0),
                          c(0,  0,  1,    0, 0, 0, 0, 0))
          ci.opt <- c(-0.99999999, 1, 0, 0)
        }
        if (K.two > 1) {
          if (weighting.two == "beta.restricted") {
            lf <- function(p) {
              llh_mf_g(df = df_llh,
                     y = ret,
                     x = covariate,
                     low.freq = low.freq,
                     mu = p["mu"],
                     omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                     alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                     m = p["m"], theta = p["theta"],
                     w1 = 1, w2 = p["w2"], g_zero = g_zero, K = K,
                     x.two = covariate.two,
                     K.two = K.two, low.freq.two = low.freq.two,
                     theta.two = p["theta.two"], w1.two = 1, w2.two = p["w2.two"])
            }
            par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                           m = 0, theta = 0, w2 = 3, theta.two = 0, w2.two = 3)
            ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0, 0),
                            c(0,  0,  0,    0, 0, 0, 1, 0, 0),
                            c(0,  1,  0,    0, 0, 0, 0, 0, 0),
                            c(0,  0,  1,    0, 0, 0, 0, 0, 0),
                            c(0,  0,  0,    0, 0, 0, 0, 0, 1))
            ci.opt <- c(-0.99999999, 1, 0, 0, 1)
          }
          if (weighting.two != "beta.restricted") {
            stop("Weighting scheme for second variable can only be beta.restricted.")
          }
        }
      }
      
      if (weighting == "beta.unrestricted" & is.null(K.two) == TRUE){
        lf <- function(p) {
          llh_mf_g(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                 m = p["m"], theta = p["theta"], w1 = p["w1"], w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                       m = 0, theta = 0, w1 = 1.0000001, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0),
                        c(0,  0,  0,  0,   0, 0, 1, 0),
                        c(0,  0,  0,  0,   0, 0, 0, 1),
                        c(0,  1,  0,  0,   0, 0, 0, 0),
                        c(0,  0,  1,  0,   0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 1, 0, 0)
      }
      
      if (weighting == "beta.unrestricted" & is.null(weighting.two) == FALSE) {
        if (weighting.two == "beta.restricted") {
          lf <- function(p) {
            llh_mf_g(df = df_llh,
                   y = ret,
                   x = covariate,
                   low.freq = low.freq,
                   mu = p["mu"],
                   omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                   alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                   m = p["m"], theta = p["theta"],
                   w1 = p["w1"], w2 = p["w2"], g_zero = g_zero, K = K,
                   x.two = covariate.two,
                   K.two = K.two, low.freq.two = low.freq.two,
                   theta.two = p["theta.two"], w1.two = 1, w2.two = p["w2.two"])
          }
          par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                         m = 0, theta = 0, w1 = 1.00000001, w2 = 3, theta.two = 0, w2.two = 3)
          ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 1, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 0, 1, 0, 0),
                          c(0,  1,  0,    0, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  1,    0, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 0, 0, 0, 1))
          ci.opt <- c(-0.99999999, 1, 1, 0, 0, 1)
        }
      }
    }
    
    if (gamma == FALSE) {
      
      if (weighting == "beta.restricted") {
        lf <- function(p) {
          llh_mf_g(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"], omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"], beta = p["beta"], gamma = 0,
                 m = p["m"], theta = p["theta"], w1 = 1, w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, 0, 0, 0),
                        c(0, 0, 0,  0, 0, 1),
                        c(0, 1, 0, 0,  0, 0),
                        c(0, 0, 1, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 0, 0)
      }
      
      if (weighting == "beta.unrestricted") {
        lf <- function(p) {
          llh_mf_g(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = 0,
                 m = p["m"],
                 theta = p["theta"],
                 w1 = p["w1"],
                 w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, w1 = 1.00000001, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, 0, 0, 0, 0),
                        c(0,  0,  0, 0, 0, 1, 0),
                        c(0,  0,  0, 0, 0, 0, 1),
                        c(0,  1,  0, 0, 0, 0, 0),
                        c(0,  0,  1, 0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 1, 0, 0)
      }
      
    }
    
    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }
    
    p.e.nlminb <- constrOptim(theta = par.start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    p.e.nlminb$value <- -p.e.nlminb$value
    
    if (multi.start == TRUE && gamma == TRUE) {
      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"] + theta["gamma"]/2 >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.two) != "try-error" && -p.e.nlminb.two$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.two
        p.e.nlminb$value <- -p.e.nlminb$value
      }
      
      par.max.lik.nr <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NR")}, silent = TRUE)
      if (class(par.max.lik.nr) != "try-error" && par.max.lik.nr$maximum > p.e.nlminb$value &&
          par.max.lik.nr$estimate["w2"] >= 1 &&
          par.max.lik.nr$estimate["alpha"] + par.max.lik.nr$estimate["beta"] + par.max.lik.nr$estimate["gamma"] / 2 < 1 &&
          par.max.lik.nr$estimate["alpha"] >= 0 && par.max.lik.nr$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nr$estimate
        p.e.nlminb$value <- par.max.lik.nr$maximum
      }
      par.max.lik.nm <- try({maxLik(logLik = function(x) -lf(x), start = par.start, method = "NM")}, silent = TRUE)
      if (class(par.max.lik.nm) != "try-error" && par.max.lik.nm$maximum > p.e.nlminb$value &&
          par.max.lik.nm$estimate["w2"] >= 1 &&
          par.max.lik.nm$estimate["alpha"] + par.max.lik.nm$estimate["beta"] + par.max.lik.nm$estimate["gamma"] / 2 < 1 &&
          par.max.lik.nm$estimate["alpha"] >= 0 && par.max.lik.nm$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nm$estimate
        p.e.nlminb$value <- par.max.lik.nm$maximum
      }
      
      p.e.nlminb.three <- try({
        suppressWarnings(optim(par = par.start, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"] + theta["gamma"]/2  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.three) != "try-error" && -p.e.nlminb.three$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.three
        p.e.nlminb$value <- -p.e.nlminb$value
      }
      
    }
    
    if (multi.start == TRUE && gamma == FALSE) {
      
      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"]  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.two) != "try-error" && -p.e.nlminb.two$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.two
        p.e.nlminb$value <- -p.e.nlminb$value
      }
      
      par.max.lik.nr <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NR")}, silent = TRUE)
      if (class(par.max.lik.nr) != "try-error" && par.max.lik.nr$maximum > p.e.nlminb$value &&
          par.max.lik.nr$estimate["w2"] >= 1 &&
          par.max.lik.nr$estimate["alpha"] + par.max.lik.nr$estimate["beta"] < 1 &&
          par.max.lik.nr$estimate["alpha"] >= 0 && par.max.lik.nr$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nr$estimate
        p.e.nlminb$value <- par.max.lik.nr$maximum
      }
      par.max.lik.nm <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NM")}, silent = TRUE)
      if (class(par.max.lik.nm) != "try-error" && par.max.lik.nm$maximum > p.e.nlminb$value &&
          par.max.lik.nm$estimate["w2"] >= 1 &&
          par.max.lik.nm$estimate["alpha"] + par.max.lik.nm$estimate["beta"] < 1 &&
          par.max.lik.nm$estimate["alpha"] >= 0 && par.max.lik.nm$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nm$estimate
        p.e.nlminb$value <- par.max.lik.nm$maximum
      }
      
      p.e.nlminb.three <- try({
        suppressWarnings(optim(par = par.start, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"]  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.three) != "try-error" && -p.e.nlminb.three$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.three
        p.e.nlminb$value <- -p.e.nlminb$value
      }
    }
    par <- p.e.nlminb$par
    
    if (weighting == "beta.restricted") {
      if (is.null(x.two) == FALSE) {
        if (K.two > 1) {
          tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                  w1 = 1, w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                  x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                  low.freq.two = low.freq.two,
                                  w1.two = 1, w2.two = par["w2.two"])$tau
        } else {
          tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                  w1 = 1, w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                  x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                  low.freq.two = low.freq.two,
                                  w1.two = 1, w2.two = 1)$tau
        }
        
      } else {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = 1, w2 = par["w2"],
                                theta = par["theta"],
                                m = par["m"], K = K)$tau
      }
      
      tau_forecast <-
        exp(sum_tau_fcts(m = par["m"],
                         i = K + 1,
                         theta = par["theta"],
                         phivar = calculate_phi(w1 = 1, w2 = par["w2"], K = K),
                         covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                         K = K))
      
      if (is.null(x.two) == FALSE) {
        if (K.two > 1) {
          tau_forecast <-
            tau_forecast *
            exp(sum_tau_fcts(m = 0,
                             i = K.two + 1,
                             theta = par["theta.two"],
                             phivar = calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two),
                             covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                             K = K.two))
        } else {
          tau_forecast <-
            tau_forecast *
            exp(sum_tau_fcts(m = 0,
                             i = K.two + 1,
                             theta = par["theta.two"],
                             phivar = calculate_phi(w1 = 1, w2 = 1, K = K.two),
                             covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                             K = K.two))
        }
        
      }
    }
    
    if (weighting == "beta.unrestricted") {
      if (is.null(x.two) == FALSE) {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = par["w1"], w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                low.freq.two = low.freq.two,
                                w1.two = 1, w2.two = par["w2.two"])$tau
      } else {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = par["w1"], w2 = par["w2"],
                                theta = par["theta"],
                                m = par["m"], K = K)$tau
      }
      
      
      tau_forecast <-
        exp(sum_tau_fcts(m = par["m"],
                         i = K + 1,
                         theta = par["theta"],
                         phivar = calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K),
                         covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                         K = K))
      if (is.null(x.two) == FALSE) {
        tau_forecast <-
          tau_forecast *
          exp(sum_tau_fcts(m = 0,
                           i = K.two + 1,
                           theta = par["theta.two"],
                           phivar = calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two),
                           covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                           K = K.two))
      }
    }
    
    returns <- unlist(data[y])
    
    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = 0,
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
    }
    
    if ((var.ratio.freq %in% c("date", low.freq)) == FALSE) {
      if (is.null(x.two) == TRUE) {
        df.fitted <- cbind(data[c("date", y, low.freq, x, var.ratio.freq)], g = g, tau = tau)
      } else {
        df.fitted <- cbind(data[c("date", y, low.freq, x, low.freq.two, x.two, var.ratio.freq)], g = g, tau = tau)
      }
      
    } else {
      if (is.null(x.two) == TRUE) {
        df.fitted <- cbind(data[c("date", y, low.freq, x)], g = g, tau = tau)
      } else {
        df.fitted <- cbind(data[c("date", y, low.freq, x, low.freq.two, x.two)], g = g, tau = tau)
      }
      
    }
    
    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
    
  }
  df.fitted$date <- as.Date(date_backup)
  # Standard errors --------------------------------------------------------------------------------
  # inv_hessian <- try({
  #   solve(-optimHess(par = par, fn = function (theta) {
  #       if( is.na(sum(lf(theta))) == TRUE) {
  #         10000000
  #       } else {
  #         sum(lf(theta))
  #       }
  #     }))
  #   }, silent = TRUE)
  
  inv_hessian <- try({
    solve(-suppressWarnings(hessian(x = par, func = function (theta) {
      if( is.na(sum(lf(theta))) == TRUE) {
        0
      } else {
        -sum(lf(theta))
      }
    })))
  }, silent = TRUE)
  
  opg.std.err <- try({sqrt(diag(solve(crossprod(jacobian(func = function(theta) -lf(theta), x = par)))))},
                     silent = TRUE)
  if (class(opg.std.err)[1] == "try-error") {
    warning("Inverting the OPG matrix failed. No OPG standard errors calculated.")
    opg.std.err <- NA
  } else {
    opg.std.err <- opg.std.err * sqrt((mean(df.fitted$residuals^4, na.rm = TRUE) - 1) / 2)
  }
  
  if (class(inv_hessian)[1] == "try-error") {
    warning("Inverting the Hessian matrix failed. No robust standard errors calculated. Possible workaround: Multiply returns by 100.")
    rob.std.err <- NA
  } else {
    rob.std.err <- sqrt(diag(inv_hessian %*% crossprod(jacobian(func = lf, x = par)) %*% inv_hessian))
  }
  
  # Output -----------------------------------------------------------------------------------------
  output <-
    list(par = par,
         std.err = rob.std.err,
         broom.mgarch = data.frame(term = names(par),
                                   estimate = par,
                                   rob.std.err = rob.std.err,
                                   p.value = 2 * (1 - pnorm(unlist(abs(par/rob.std.err)))),
                                   opg.std.err = opg.std.err,
                                   opg.p.value = 2 * (1 - pnorm(unlist(abs(par/opg.std.err))))),
         tau = tau,
         g = g,
         df.fitted = df.fitted,
         K = K,
         weighting.scheme = weighting,
         llh = p.e.nlminb$value,
         bic = log(sum(!is.na(tau))) * length(par) - 2 * (p.e.nlminb$value),
         y = y,
         optim = p.e.nlminb)
  
  if (is.null(x.two) == FALSE) {
    output$K.two <- K.two
    output$weighting.scheme.two <- weighting.two
  }
  if (K == 0) {
    output$tau.forecast <- exp(par["m"])
  }
  
  
  # Additional output if there is a long-term component (K > 0) -------------------------------------
  if (K > 0) {
    output$variance.ratio <- 100 *
      var(log(aggregate(df.fitted$tau, by = df.fitted[var.ratio.freq],
                        FUN = mean)[,2]),
          na.rm = TRUE) /
      var(log(aggregate(df.fitted$tau * df.fitted$g, by = df.fitted[var.ratio.freq],
                        FUN = mean)[,2]),
          na.rm = TRUE)
    output$tau.forecast <- tau_forecast
    
    if (weighting == "beta.restricted") {
      output$est.weighting <- calculate_phi(1, w2 = par["w2"], K = K)
    }
    if (weighting == "beta.unrestricted") {
      output$est.weighting <- calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K)
    }
    if (is.null(x.two) == FALSE) {
      if (K.two > 1) {
        output$est.weighting.two <- calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two)
      }
    }
    
  }
  
  # Add class mfGARCH for employing generic functions
  class(output) <- "mfGARCH"
  output
}

# esempio per vedere se funziona con la distrubuzione Normale
v<-data.frame(df_mfgarch[,1:11])
v<-v[1:11644,]
fit_g_dhousing<-fit_mfgarch_g(data = v, y = "return", x = "dhousing", K = 36, low.freq = "year_month")


### CASO t-Student ###
llh_mf_t <-
  function(df, x, y, low.freq, mu, omega, alpha, beta, gamma,
           m, theta, w1 = 1, w2 = 1, g_zero, K,
           x.two = NULL, K.two = NULL, theta.two = NULL,
           low.freq.two = NULL, w1.two = NULL, w2.two = NULL) {
    
    if (is.null(x.two) == FALSE) {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K,
                              x.two = x.two, K.two = K.two, theta.two = theta.two,
                              low.freq.two = low.freq.two,
                              w1.two = w1.two, w2.two = w2.two)$tau
    } else {
      tau <- calculate_tau_mf(df = df, x = x, low.freq = low.freq,
                              w1 = w1, w2 = w2, theta = theta, m = m, K = K)$tau
    }
    
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
    ll<-ll <- -dstd(ret, mean = mu, sd = g*tau, nu = 30, log = TRUE)
    }
  }

llh_simple_t <- function(y, mu, alpha, beta, gamma, m, g_zero) {
  omega <- 1 - alpha - beta - gamma / 2
  ret <- y
  ret_std <- (ret - mu)/sqrt(exp(m))
  g <- calculate_g(omega = omega, alpha = alpha, beta = beta, gamma = gamma,
                   returns = ret_std, g0 = g_zero)
  ll<-ll <- -dstd(ret, mean = mu, sd = g*exp(m), nu = 30, log = TRUE)
}

# funzione generica che mi stima parametri assumendo t-Student
fit_mfgarch_t <- function(data, y, x = NULL, K = NULL, low.freq = "date", var.ratio.freq = NULL, gamma = TRUE, weighting = "beta.restricted", x.two = NULL, K.two = NULL, low.freq.two = NULL, weighting.two = NULL, multi.start = FALSE, control = list(par.start = NULL)) {
  
  print("For ensuring numerical stability of the parameter optimization and inversion of the Hessian, it is best to multiply log returns by 100.")
  
  if (is.null(weighting.two) == FALSE) {
    if (weighting.two != "beta.restricted") {
      stop("Right now, only beta.restricted weighting scheme is employed for the second covariate.")
    }
  }
  
  if (is.null(x.two) == FALSE) {
    weighting.two <- "beta.restricted"
  }
  
  if (is.null(x.two) == FALSE && gamma == FALSE) {
    stop("Regarding two covariates, only asymmetric GJR-GARCH component is implemented.")
  }
  
  # if (K == 1 && is.null(x.two) == FALSE && K.two != 1) {
  #   stop("Regarding two covariates, only K.two = 1 is implemented if K = 1.")
  # }
  if (is.null(K.two) == FALSE) {
    if (K == 1 & K.two > 1) {
      stop("Regarding two covariates with one of them being equal to one, only K.two = 1 is implemented.")
    }
  }
  
  if (is.null(x.two) == FALSE) {
    print("Specifying two covariates may lead to long estimation times.")
  }
  
  if (weighting %in% c("beta.restricted", "beta.unrestricted") == FALSE) {
    stop("Incorrect weighting scheme specified - options are \"beta.restricted\" and \"beta.unrestricted\".")
  }
  if (gamma %in% c(TRUE, FALSE) == FALSE) {
    stop("Gamma can't be anything different than TRUE or FALSE.")
  }
  if ("date" %in% colnames(data) == FALSE) {
    stop("No date column.")
  }
  if (inherits(data$date, 'Date') == FALSE) {
    stop(paste0("Supplied date column is not of format 'Date'. It is of class '", class(data$date), "'."))
  }
  if (inherits(data[[low.freq]], 'Date') == FALSE) {
    stop(paste0("Supplied low.freq column is not of format 'Date'. It is of class '", class(data[[low.freq]]), "'."))
  }
  if (is.null(x) == FALSE && K == 0) {
    warning("You specified an external covariate x but chose K = 0 - simple GARCH is estimated (K = 0).")
  }
  
  if (is.null(x) == TRUE) {
    warning("No external covariate x is specified - simple GARCH is estimated (K = 0).")
    x <- "date"
    K <- 0
  }
  if (is.null(K) == TRUE) {
    warning("No K is specified - simple GARCH is estimated (K = 0).")
    x <- "date"
    K <- 0
  }
  if (K < 0 || K %% 1 != 0) {
    stop("K can't be smaller than 0 and has to be an integer.")
  }
  if (dim(unique(data[c(x, low.freq)]))[1] > dim(unique(data[c(low.freq)]))[1]) {
    stop("There is more than one unique observation per low frequency entry.")
  }
  # if ((is.null(x) == TRUE && (is.null(K) == TRUE)) || K == 0) {
  #   K <- 0
  # }
  if (y %in% colnames(data) == FALSE) {
    stop(paste("There is no variable in your data frame with name ", y, "."))
  }
  if (x %in% colnames(data) == FALSE && is.null(x) != FALSE) {
    stop(paste("There is no variable in your data frame with name ", x, "."))
  }
  if (low.freq %in% colnames(data) == FALSE) {
    stop(paste("There is no low freq. variable in your data frame with name ", low.freq, "."))
  }
  if ("tau" %in% colnames(data) == TRUE) {
    stop("There may not be a column named tau - it will be part of df.fitted")
  }
  if ("g" %in% colnames(data) == TRUE) {
    stop("There may not be a column named g - it will be part of df.fitted")
  }
  if (is.null(x) == TRUE) {
    if (sum(is.na(data[[y]]) == TRUE) > 0) {
      stop(paste0("Column ", y, "contains NAs."))
    }
  } else {
    if (sum(is.na(data[[y]]) == TRUE) > 0 | sum(is.na(data[[x]]) == TRUE) > 0) {
      stop(paste0("Either column ", y, " or column ", x, "includes NAs."))
    }
  }
  if (length(unlist(unique(data[["date"]]))) != dim(data)[1]) {
    stop("There is more than one observation per high frequency (presumably date).")
  }
  if (is.null(var.ratio.freq) == FALSE) {
    if (var.ratio.freq %in% colnames(data) == FALSE) {
      stop(paste0("There is no var.ratio.freq column with name ", var.ratio.freq, "."))
    }
  }
  
  # Order by high frequency variable
  data <- data[order(data$date), ]
  # Deprecated dplyr version
  #data <- dplyr::arrange_(data, "date")
  # We store date in new variable because computation on integerized date seemed to be faster
  date_backup <- data[["date"]]
  data["date"] <- as.numeric(unlist(data["date"]))
  
  if (is.null(var.ratio.freq) == TRUE) {
    var.ratio.freq <- low.freq
    print(paste0("No frequency specified for calculating the variance ratio - default: low.freq = ", low.freq))
  }
  
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
  
  g_zero <- var(unlist(data[[y]]))
  ret <- data[[y]]
  
  # Parameter estimation
  if (K == 0) {
    if (gamma == TRUE) {
      lf <- function(p) {
        llh_simple_t(y = ret,
                   mu = p["mu"],
                   alpha = p["alpha"],
                   beta = p["beta"],
                   gamma = p["gamma"],
                   m = p["m"],
                   g_zero = g_zero)
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0)
      ui.opt <- rbind(c(0, -1, -1, -1/2, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
      ci.opt <- c(-0.99999, 0, 0)
    } else {
      lf <- function(p) {
        llh_simple_t(y = ret,
                   mu = p["mu"],
                   alpha = p["alpha"],
                   beta = p["beta"],
                   gamma = 0,
                   m = p["m"],
                   g_zero = g_zero)
        
      }
      par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0)
      ui.opt <- rbind(c(0, -1, -1, 0), c(0, 1, 0, 0), c(0, 0, 1, 0))
      ci.opt <- c(-0.99999, 0, 0)
    }
    
    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }
    
    p.e.nlminb <- constrOptim(theta = par.start,
                              f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    
    if (multi.start == TRUE && gamma == TRUE) {
      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.two) == "try-error") {
        print("Second-step BFGS optimization failed. Fallback: First-stage Nelder-Mead estimate.")
      } else {
        if (p.e.nlminb.two$value < p.e.nlminb$value) {
          p.e.nlminb <- p.e.nlminb.two
        }
      }
    }
    
    p.e.nlminb$value <- -p.e.nlminb$value
    
    par <- p.e.nlminb$par
    returns <- as.numeric(unlist(data[[y]]))
    tau <- rep(exp(par["m"]), times = length(returns))
    
    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] -  par["gamma"]/2,
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
      tau <- rep(exp(par["m"]), times = length(g))
      
    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = 0,
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
      tau <- rep(exp(par["m"]), times = length(g))
    }
    
    if ((var.ratio.freq %in% c("date", "low.freq")) == FALSE) {
      df.fitted <- cbind(data[c("date", y, var.ratio.freq)], g = g, tau = tau)
    } else {
      df.fitted <- cbind(data[c("date", y)], g = g, tau = tau)
    }
    
    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
  } else { # if K > 0 we get the covariate series
    covariate <- unlist(unique(data[c(low.freq, x)])[x])
    
    if (is.null(x.two) == FALSE) {
      covariate.two <- unlist(unique(data[c(low.freq.two, x.two)])[x.two])
    }
  }
  
  if (K == 1) {
    if (is.null(K.two) == FALSE) {
      if (gamma == TRUE) {
        lf <- function(p) {
          llh_mf_t(df = df_llh, y = ret, x = covariate, low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = p["gamma"],
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K,
                 x.two = covariate.two,
                 K.two = K.two, low.freq.two = low.freq.two,
                 theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0, theta.two = 0)
        ui_opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0), c(0, 1, 0, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      } else {
        lf <- function(p) {
          llh_mf_t(df = df_llh,
                 y = ret, x = covariate, low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = 0,
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K,
                 x.two = covariate.two,
                 K.two = K.two, low.freq.two = low.freq.two,
                 theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, theta.two = 0)
        ui_opt <- rbind(c(0, -1, -1,  0, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      }
      
    } else {
      
      if (gamma == TRUE) {
        lf <- function(p) {
          llh_mf_t(df = df_llh, y = ret, x = covariate, low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = p["gamma"],
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04, m = 0, theta = 0)
        ui_opt <- rbind(c(0, -1, -1, -1/2, 0, 0), c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      } else {
        lf <- function(p) {
          llh_mf_t(df = df_llh,
                 y = ret, x = covariate, low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = 0,
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1, w2 = 1, g_zero = g_zero, K = K)
        }
        par_start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0)
        ui_opt <- rbind(c(0, -1, -1,  0, 0), c(0, 1, 0, 0, 0), c(0, 0, 1, 0, 0))
        ci_opt <- c(-0.99999, 0, 0)
      }
    }
    
    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }
    
    p.e.nlminb <- constrOptim(theta = par_start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui_opt, ci = ci_opt, hessian = FALSE)
    par <- p.e.nlminb$par
    p.e.nlminb$value <- -p.e.nlminb$value
    
    if (is.null(x.two) == FALSE) {
      tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                              w1 = 1, w2 = 1, theta = par["theta"], m = par["m"], K = K,
                              x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                              low.freq.two = low.freq.two,
                              w1.two = 1, w2.two = 1)$tau
    } else {
      tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                              theta = par["theta"], m = par["m"], w1 = 1, w2 = 1, K = K)$tau
    }
    
    tau_forecast <-
      exp(sum_tau_fcts(m = par["m"],
                       i = K + 1,
                       theta = par["theta"],
                       phivar = calculate_phi(w1 = 1, w2 = 1, K = K),
                       covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                       K = K))
    
    if (is.null(x.two) == FALSE) {
      tau_forecast <-
        tau_forecast *
        exp(sum_tau_fcts(m = 0,
                         i = K.two + 1,
                         theta = par["theta.two"],
                         phivar = calculate_phi(w1 = 1, w2 = 1, K = K.two),
                         covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                         K = K.two))
    }
    
    
    returns <- unlist(data[y])
    
    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                         alpha = par["alpha"], beta = par["beta"], gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))), g0 = g_zero))
    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],
                         alpha = par["alpha"], beta = par["beta"], gamma = 0,
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))), g0 = g_zero))
    }
    
    if ((var.ratio.freq %in% c("date", "low.freq")) == FALSE) {
      df.fitted <- cbind(data[c("date", y, low.freq, x, var.ratio.freq)], g = g, tau = tau)
    } else {
      df.fitted <- cbind(data[c("date", y, low.freq, x)], g = g, tau = tau)
    }
    
    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
    
  }
  
  if (K > 1) {
    if (gamma == TRUE) {
      if (weighting == "beta.restricted" & is.null(K.two) == TRUE) {
        lf <- function(p) {
          llh_mf_t(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = p["gamma"],
                 m = p["m"],
                 theta = p["theta"],
                 w1 = 1,
                 w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                       m = 0, theta = 0, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0),
                        c(0,  0,  0,    0, 0, 0, 1),
                        c(0,  1,  0,    0, 0, 0, 0),
                        c(0,  0,  1,    0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 0, 0)
      }
      if (weighting == "beta.restricted" & is.null(K.two) == FALSE) {
        if (K.two == 1) {
          lf <- function(p) {
            llh_mf_t(df = df_llh,
                   y = ret,
                   x = covariate,
                   low.freq = low.freq,
                   mu = p["mu"],
                   omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                   alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                   m = p["m"], theta = p["theta"],
                   w1 = 1, w2 = p["w2"], g_zero = g_zero, K = K,
                   x.two = covariate.two,
                   K.two = 1, low.freq.two = low.freq.two,
                   theta.two = p["theta.two"], w1.two = 1, w2.two = 1)
          }
          par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                         m = 0, theta = 0, w2 = 3, theta.two = 0)
          ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 1, 0),
                          c(0,  1,  0,    0, 0, 0, 0, 0),
                          c(0,  0,  1,    0, 0, 0, 0, 0))
          ci.opt <- c(-0.99999999, 1, 0, 0)
        }
        if (K.two > 1) {
          if (weighting.two == "beta.restricted") {
            lf <- function(p) {
              llh_mf_t(df = df_llh,
                     y = ret,
                     x = covariate,
                     low.freq = low.freq,
                     mu = p["mu"],
                     omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                     alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                     m = p["m"], theta = p["theta"],
                     w1 = 1, w2 = p["w2"], g_zero = g_zero, K = K,
                     x.two = covariate.two,
                     K.two = K.two, low.freq.two = low.freq.two,
                     theta.two = p["theta.two"], w1.two = 1, w2.two = p["w2.two"])
            }
            par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                           m = 0, theta = 0, w2 = 3, theta.two = 0, w2.two = 3)
            ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0, 0),
                            c(0,  0,  0,    0, 0, 0, 1, 0, 0),
                            c(0,  1,  0,    0, 0, 0, 0, 0, 0),
                            c(0,  0,  1,    0, 0, 0, 0, 0, 0),
                            c(0,  0,  0,    0, 0, 0, 0, 0, 1))
            ci.opt <- c(-0.99999999, 1, 0, 0, 1)
          }
          if (weighting.two != "beta.restricted") {
            stop("Weighting scheme for second variable can only be beta.restricted.")
          }
        }
      }
      
      if (weighting == "beta.unrestricted" & is.null(K.two) == TRUE){
        lf <- function(p) {
          llh_mf_t(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                 alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                 m = p["m"], theta = p["theta"], w1 = p["w1"], w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                       m = 0, theta = 0, w1 = 1.0000001, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0),
                        c(0,  0,  0,  0,   0, 0, 1, 0),
                        c(0,  0,  0,  0,   0, 0, 0, 1),
                        c(0,  1,  0,  0,   0, 0, 0, 0),
                        c(0,  0,  1,  0,   0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 1, 0, 0)
      }
      
      if (weighting == "beta.unrestricted" & is.null(weighting.two) == FALSE) {
        if (weighting.two == "beta.restricted") {
          lf <- function(p) {
            llh_mf_t(df = df_llh,
                   y = ret,
                   x = covariate,
                   low.freq = low.freq,
                   mu = p["mu"],
                   omega = 1 - p["alpha"] - p["beta"] - p["gamma"]/2,
                   alpha = p["alpha"], beta = p["beta"], gamma = p["gamma"],
                   m = p["m"], theta = p["theta"],
                   w1 = p["w1"], w2 = p["w2"], g_zero = g_zero, K = K,
                   x.two = covariate.two,
                   K.two = K.two, low.freq.two = low.freq.two,
                   theta.two = p["theta.two"], w1.two = 1, w2.two = p["w2.two"])
          }
          par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, gamma = 0.04,
                         m = 0, theta = 0, w1 = 1.00000001, w2 = 3, theta.two = 0, w2.two = 3)
          ui.opt <- rbind(c(0, -1, -1, -1/2, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 1, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 0, 1, 0, 0),
                          c(0,  1,  0,    0, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  1,    0, 0, 0, 0, 0, 0, 0),
                          c(0,  0,  0,    0, 0, 0, 0, 0, 0, 1))
          ci.opt <- c(-0.99999999, 1, 1, 0, 0, 1)
        }
      }
    }
    
    if (gamma == FALSE) {
      
      if (weighting == "beta.restricted") {
        lf <- function(p) {
          llh_mf_t(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"], omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"], beta = p["beta"], gamma = 0,
                 m = p["m"], theta = p["theta"], w1 = 1, w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, 0, 0, 0),
                        c(0, 0, 0,  0, 0, 1),
                        c(0, 1, 0, 0,  0, 0),
                        c(0, 0, 1, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 0, 0)
      }
      
      if (weighting == "beta.unrestricted") {
        lf <- function(p) {
          llh_mf_t(df = df_llh,
                 y = ret,
                 x = covariate,
                 low.freq = low.freq,
                 mu = p["mu"],
                 omega = 1 - p["alpha"] - p["beta"],
                 alpha = p["alpha"],
                 beta = p["beta"],
                 gamma = 0,
                 m = p["m"],
                 theta = p["theta"],
                 w1 = p["w1"],
                 w2 = p["w2"],
                 g_zero = g_zero,
                 K = K)
        }
        par.start <- c(mu = 0, alpha = 0.02, beta = 0.85, m = 0, theta = 0, w1 = 1.00000001, w2 = 3)
        ui.opt <- rbind(c(0, -1, -1, 0, 0, 0, 0),
                        c(0,  0,  0, 0, 0, 1, 0),
                        c(0,  0,  0, 0, 0, 0, 1),
                        c(0,  1,  0, 0, 0, 0, 0),
                        c(0,  0,  1, 0, 0, 0, 0))
        ci.opt <- c(-0.99999999, 1, 1, 0, 0)
      }
      
    }
    
    if(is.null(control$par.start) == FALSE) {
      par.start <- control$par.start
    }
    
    p.e.nlminb <- constrOptim(theta = par.start, f = function(theta) { sum(lf(theta)) },
                              grad = NULL, ui = ui.opt, ci = ci.opt, hessian = FALSE)
    p.e.nlminb$value <- -p.e.nlminb$value
    
    if (multi.start == TRUE && gamma == TRUE) {
      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"] + theta["gamma"]/2 >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.two) != "try-error" && -p.e.nlminb.two$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.two
        p.e.nlminb$value <- -p.e.nlminb$value
      }
      
      par.max.lik.nr <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NR")}, silent = TRUE)
      if (class(par.max.lik.nr) != "try-error" && par.max.lik.nr$maximum > p.e.nlminb$value &&
          par.max.lik.nr$estimate["w2"] >= 1 &&
          par.max.lik.nr$estimate["alpha"] + par.max.lik.nr$estimate["beta"] + par.max.lik.nr$estimate["gamma"] / 2 < 1 &&
          par.max.lik.nr$estimate["alpha"] >= 0 && par.max.lik.nr$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nr$estimate
        p.e.nlminb$value <- par.max.lik.nr$maximum
      }
      par.max.lik.nm <- try({maxLik(logLik = function(x) -lf(x), start = par.start, method = "NM")}, silent = TRUE)
      if (class(par.max.lik.nm) != "try-error" && par.max.lik.nm$maximum > p.e.nlminb$value &&
          par.max.lik.nm$estimate["w2"] >= 1 &&
          par.max.lik.nm$estimate["alpha"] + par.max.lik.nm$estimate["beta"] + par.max.lik.nm$estimate["gamma"] / 2 < 1 &&
          par.max.lik.nm$estimate["alpha"] >= 0 && par.max.lik.nm$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nm$estimate
        p.e.nlminb$value <- par.max.lik.nm$maximum
      }
      
      p.e.nlminb.three <- try({
        suppressWarnings(optim(par = par.start, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"] + theta["gamma"]/2  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.three) != "try-error" && -p.e.nlminb.three$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.three
        p.e.nlminb$value <- -p.e.nlminb$value
      }
      
    }
    
    if (multi.start == TRUE && gamma == FALSE) {
      
      p.e.nlminb.two <- try({
        suppressWarnings(optim(par = p.e.nlminb$par, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"]  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.two) != "try-error" && -p.e.nlminb.two$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.two
        p.e.nlminb$value <- -p.e.nlminb$value
      }
      
      par.max.lik.nr <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NR")}, silent = TRUE)
      if (class(par.max.lik.nr) != "try-error" && par.max.lik.nr$maximum > p.e.nlminb$value &&
          par.max.lik.nr$estimate["w2"] >= 1 &&
          par.max.lik.nr$estimate["alpha"] + par.max.lik.nr$estimate["beta"] < 1 &&
          par.max.lik.nr$estimate["alpha"] >= 0 && par.max.lik.nr$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nr$estimate
        p.e.nlminb$value <- par.max.lik.nr$maximum
      }
      par.max.lik.nm <- try({maxLik(logLik = function(x) - lf(x), start = par.start, method = "NM")}, silent = TRUE)
      if (class(par.max.lik.nm) != "try-error" && par.max.lik.nm$maximum > p.e.nlminb$value &&
          par.max.lik.nm$estimate["w2"] >= 1 &&
          par.max.lik.nm$estimate["alpha"] + par.max.lik.nm$estimate["beta"] < 1 &&
          par.max.lik.nm$estimate["alpha"] >= 0 && par.max.lik.nm$estimate["beta"] >= 0) {
        p.e.nlminb$par <- par.max.lik.nm$estimate
        p.e.nlminb$value <- par.max.lik.nm$maximum
      }
      
      p.e.nlminb.three <- try({
        suppressWarnings(optim(par = par.start, fn = function (theta) {
          if( is.na(sum(lf(theta))) == TRUE | theta["alpha"] < 0 | theta["alpha"] + theta["beta"]  >= 1 | theta["w2"] < 1) {
            NA
          } else {
            sum(lf(theta))
          }
        }, method = "BFGS"))}, silent = TRUE)
      
      if (class(p.e.nlminb.three) != "try-error" && -p.e.nlminb.three$value > p.e.nlminb$value) {
        p.e.nlminb <- p.e.nlminb.three
        p.e.nlminb$value <- -p.e.nlminb$value
      }
    }
    par <- p.e.nlminb$par
    
    if (weighting == "beta.restricted") {
      if (is.null(x.two) == FALSE) {
        if (K.two > 1) {
          tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                  w1 = 1, w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                  x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                  low.freq.two = low.freq.two,
                                  w1.two = 1, w2.two = par["w2.two"])$tau
        } else {
          tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                  w1 = 1, w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                  x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                  low.freq.two = low.freq.two,
                                  w1.two = 1, w2.two = 1)$tau
        }
        
      } else {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = 1, w2 = par["w2"],
                                theta = par["theta"],
                                m = par["m"], K = K)$tau
      }
      
      tau_forecast <-
        exp(sum_tau_fcts(m = par["m"],
                         i = K + 1,
                         theta = par["theta"],
                         phivar = calculate_phi(w1 = 1, w2 = par["w2"], K = K),
                         covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                         K = K))
      
      if (is.null(x.two) == FALSE) {
        if (K.two > 1) {
          tau_forecast <-
            tau_forecast *
            exp(sum_tau_fcts(m = 0,
                             i = K.two + 1,
                             theta = par["theta.two"],
                             phivar = calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two),
                             covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                             K = K.two))
        } else {
          tau_forecast <-
            tau_forecast *
            exp(sum_tau_fcts(m = 0,
                             i = K.two + 1,
                             theta = par["theta.two"],
                             phivar = calculate_phi(w1 = 1, w2 = 1, K = K.two),
                             covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                             K = K.two))
        }
        
      }
    }
    
    if (weighting == "beta.unrestricted") {
      if (is.null(x.two) == FALSE) {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = par["w1"], w2 = par["w2"], theta = par["theta"], m = par["m"], K = K,
                                x.two = covariate.two, K.two = K.two, theta.two = par["theta.two"],
                                low.freq.two = low.freq.two,
                                w1.two = 1, w2.two = par["w2.two"])$tau
      } else {
        tau <- calculate_tau_mf(df = data, x = covariate, low.freq = low.freq,
                                w1 = par["w1"], w2 = par["w2"],
                                theta = par["theta"],
                                m = par["m"], K = K)$tau
      }
      
      
      tau_forecast <-
        exp(sum_tau_fcts(m = par["m"],
                         i = K + 1,
                         theta = par["theta"],
                         phivar = calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K),
                         covariate = c(tail(unlist(unique(data[c(x, low.freq)])[x]), K), NA),
                         K = K))
      if (is.null(x.two) == FALSE) {
        tau_forecast <-
          tau_forecast *
          exp(sum_tau_fcts(m = 0,
                           i = K.two + 1,
                           theta = par["theta.two"],
                           phivar = calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two),
                           covariate = c(tail(unlist(unique(data[c(x.two, low.freq.two)])[x.two]), K.two), NA),
                           K = K.two))
      }
    }
    
    returns <- unlist(data[y])
    
    if (gamma == TRUE) {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"] - par["gamma"]/2,
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = par["gamma"],
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
    } else {
      g <- c(rep(NA, times = sum(is.na((returns - par["mu"])/sqrt(tau)))),
             calculate_g(omega = 1 - par["alpha"] - par["beta"],
                         alpha = par["alpha"],
                         beta = par["beta"],
                         gamma = 0,
                         as.numeric(na.exclude((returns - par["mu"])/sqrt(tau))),
                         g0 = g_zero))
    }
    
    if ((var.ratio.freq %in% c("date", low.freq)) == FALSE) {
      if (is.null(x.two) == TRUE) {
        df.fitted <- cbind(data[c("date", y, low.freq, x, var.ratio.freq)], g = g, tau = tau)
      } else {
        df.fitted <- cbind(data[c("date", y, low.freq, x, low.freq.two, x.two, var.ratio.freq)], g = g, tau = tau)
      }
      
    } else {
      if (is.null(x.two) == TRUE) {
        df.fitted <- cbind(data[c("date", y, low.freq, x)], g = g, tau = tau)
      } else {
        df.fitted <- cbind(data[c("date", y, low.freq, x, low.freq.two, x.two)], g = g, tau = tau)
      }
      
    }
    
    df.fitted$residuals <- unlist((df.fitted[y] - par["mu"]) / sqrt(df.fitted$g * df.fitted$tau))
    
  }
  df.fitted$date <- as.Date(date_backup)
  # Standard errors --------------------------------------------------------------------------------
  # inv_hessian <- try({
  #   solve(-optimHess(par = par, fn = function (theta) {
  #       if( is.na(sum(lf(theta))) == TRUE) {
  #         10000000
  #       } else {
  #         sum(lf(theta))
  #       }
  #     }))
  #   }, silent = TRUE)
  
  inv_hessian <- try({
    solve(-suppressWarnings(hessian(x = par, func = function (theta) {
      if( is.na(sum(lf(theta))) == TRUE) {
        0
      } else {
        -sum(lf(theta))
      }
    })))
  }, silent = TRUE)
  
  opg.std.err <- try({sqrt(diag(solve(crossprod(jacobian(func = function(theta) -lf(theta), x = par)))))},
                     silent = TRUE)
  if (class(opg.std.err)[1] == "try-error") {
    warning("Inverting the OPG matrix failed. No OPG standard errors calculated.")
    opg.std.err <- NA
  } else {
    opg.std.err <- opg.std.err * sqrt((mean(df.fitted$residuals^4, na.rm = TRUE) - 1) / 2)
  }
  
  if (class(inv_hessian)[1] == "try-error") {
    warning("Inverting the Hessian matrix failed. No robust standard errors calculated. Possible workaround: Multiply returns by 100.")
    rob.std.err <- NA
  } else {
    rob.std.err <- sqrt(diag(inv_hessian %*% crossprod(jacobian(func = lf, x = par)) %*% inv_hessian))
  }
  
  # Output -----------------------------------------------------------------------------------------
  output <-
    list(par = par,
         std.err = rob.std.err,
         broom.mgarch = data.frame(term = names(par),
                                   estimate = par,
                                   rob.std.err = rob.std.err,
                                   p.value = 2 * (1 - pnorm(unlist(abs(par/rob.std.err)))),
                                   opg.std.err = opg.std.err,
                                   opg.p.value = 2 * (1 - pnorm(unlist(abs(par/opg.std.err))))),
         tau = tau,
         g = g,
         df.fitted = df.fitted,
         K = K,
         weighting.scheme = weighting,
         llh = p.e.nlminb$value,
         bic = log(sum(!is.na(tau))) * length(par) - 2 * (p.e.nlminb$value),
         y = y,
         optim = p.e.nlminb)
  
  if (is.null(x.two) == FALSE) {
    output$K.two <- K.two
    output$weighting.scheme.two <- weighting.two
  }
  if (K == 0) {
    output$tau.forecast <- exp(par["m"])
  }
  
  
  # Additional output if there is a long-term component (K > 0) -------------------------------------
  if (K > 0) {
    output$variance.ratio <- 100 *
      var(log(aggregate(df.fitted$tau, by = df.fitted[var.ratio.freq],
                        FUN = mean)[,2]),
          na.rm = TRUE) /
      var(log(aggregate(df.fitted$tau * df.fitted$g, by = df.fitted[var.ratio.freq],
                        FUN = mean)[,2]),
          na.rm = TRUE)
    output$tau.forecast <- tau_forecast
    
    if (weighting == "beta.restricted") {
      output$est.weighting <- calculate_phi(1, w2 = par["w2"], K = K)
    }
    if (weighting == "beta.unrestricted") {
      output$est.weighting <- calculate_phi(w1 = par["w1"], w2 = par["w2"], K = K)
    }
    if (is.null(x.two) == FALSE) {
      if (K.two > 1) {
        output$est.weighting.two <- calculate_phi(w1 = 1, w2 = par["w2.two"], K = K.two)
      }
    }
    
  }
  
  # Add class mfGARCH for employing generic functions
  class(output) <- "mfGARCH"
  output
}

fit_t_dhousing<-fit_mfgarch_t(data = v, y = "return", x = "dhousing", K = 36, low.freq = "year_month")
# commento: non vengono uguali a quelli della normale...
# non calcola l'Hessiana
fit_t_dhousing$par
fit_t_dhousing$bic

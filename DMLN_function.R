DMLN <- function(X,Y,T = 1000){
  N <- dim(X)[1]
  R <- dim(X)[2]
  P <- dim(Y)[2]
  # Parameters
  set.seed(42)
  tau <- 1
  a_p <- 0.4
  b_p <- 1.6
  a_beta <- 2
  b_beta <- 10
  # beta_0_mean <- 1
  
  # Initialize
  p_delta <- a_p / (a_p + b_p)
  # Delta_temp <- Delta
  Delta_temp <- matrix(rbern(R*P, p_delta), nrow = R, ncol = P)
  
  Beta_temp <- matrix(0, nrow = R, ncol = P)
  for (i in 1:R){
    for (j in 1:P){
      if (Delta_temp[i,j] == 1){
        Beta_temp[i,j] <- rnorm(1, 0, tau)
      }
    }
  }
  # Beta_temp <- Beta
  
  # Beta_0_temp <- Beta_0
  Beta_0_temp <- rnorm(P,0,1)
  X_design <- cbind(1, X)
  log_alpha <- X_design %*% rbind(Beta_0_temp, Beta_temp)
  Alpha <- exp(log_alpha)
  
  # Create space to store the result
  Delta_store <- array(data = NA, dim = c(T, R, P))
  Beta_store <- array(data = NA, dim = c(T, R, P))
  Beta_0_store <- matrix(NA, nrow = T, ncol = P)
  
  count <- 0
  ########## Start the MCMC Samplying ##########
  for (t in 1:T){
    if(t*100/T == count)
    {
      spintf(" %d % has been done!", count)
      count = count + 10;
    }
    for (j in 1:P){
      # Stochastic method
      # covar <- sample(seq(R),2)
      # for (r in covar){
      for (r in 1:R){
        # print(Delta_temp)
        Delta_temp[r,j] <- 1 - Delta_temp[r,j]
        if (Delta_temp[r,j] == 1){
          # add
          beta_star <- rnorm(1, 0, 1)
          # beta_star <- rtruncnorm(1,mean = 0, sd = 1, low = -5, high = 5)
          # method 1------
          # Beta_star <- Beta_temp
          # Beta_star[r,j] <- beta_star
          # Alpha_star <- exp(X_design %*% rbind(Beta_0_temp, Beta_star))
          # method 2------
          # b <- Beta_temp[,j]
          # b[r] <- beta_star
          # a <- X_design %*% c(Beta_0_temp[j],b)
          # Alpha_star <- Alpha
          # Alpha_star[,j] <- exp(a)
          # method 3------
          log_alpha_star <- log_alpha
          for (i in 1:N){
            log_alpha_star[i,j] <- log_alpha_star[i,j] + beta_star*X[i,r]
          }
          Alpha_star <- exp(log_alpha_star)
          
          r_mh <- 0
          for (l in 1:N){
            alpha_old <- Alpha[l,]
            alpha_star <- Alpha_star[l,]
            y <- Y[l,]
            r_mh <- r_mh + lgamma(sum(alpha_star))- lgamma(sum(alpha_star) + sum(y))- lgamma(sum(alpha_old))+ lgamma(sum(alpha_old) + sum(y))
            r_mh <- r_mh + lgamma(alpha_star[j] + y[j]) - lgamma(alpha_star[j]) - lgamma(alpha_old[j] + y[j]) + lgamma(alpha_old[j])
            # for (k in 1:P){
            #   r_mh <- r_mh + lgamma(alpha_star[k] + y[k]) - lgamma(alpha_star[k]) - lgamma(alpha_old[k] + y[k]) + lgamma(alpha_old[k])
            # }
          }
          # r_mh <- -0.5*log(2*pi) + a_beta*log(b_beta) - lgamma(a_beta) + lgamma(a_beta + 0.5) - (a_beta + 0.5)*log(b_beta + beta_star^2/2)
          sum_delta <- sum(Delta_temp[,j])
          # r_mh <- r_mh + lgamma(a_p + sum_delta) + lgamma(R-sum_delta + b_p) - lgamma(a_p + sum_delta - 1) - lgamma(R - sum_delta + 1 + b_p)
          # r_mh <- r_mh + lgamma(a_p+1) + lgamma(b_p)-lgamma(a_p)-lgamma(b_p+1)
          lb_star <- lbeta(a_p + sum_delta, R-sum_delta + b_p) - lbeta(a_p, b_p)
          lb_old <- lbeta(a_p + sum_delta - 1, R - sum_delta + 1 + b_p) - lbeta(a_p, b_p)
          r_mh <- r_mh + lb_star - log(1-exp(lb_old))
          r_mh <- r_mh + dnorm(beta_star, 0, 10, log = TRUE)
          r_mh <- r_mh - dnorm(beta_star,0, 1, log = TRUE)
          
          if (r_mh > log(runif(1))){
            Beta_temp[r,j] <- beta_star
            log_alpha <- log_alpha_star
            Alpha <- Alpha_star
          }else{
            Delta_temp[r,j] <- 0
          }
        }else{
          # delete
          beta_old <- Beta_temp[r,j]
          # Beta_star <- Beta_temp
          # Beta_star[r,j] <- 0
          # Alpha_star <- exp(X_design %*% rbind(Beta_0_temp, Beta_star))
          # b <- Beta_temp[,j]
          # b[r] <- 0
          # a <- X_design %*% c(Beta_0_temp[j],b)
          # Alpha_star <- Alpha
          # Alpha_star[,j] <- exp(a)
          
          log_alpha_star <- log_alpha
          for (i in 1:N){
            log_alpha_star[i,j] <- log_alpha_star[i,j] - beta_old*X[i,r]
          }
          Alpha_star <- exp(log_alpha_star)
          
          r_mh <- 0
          for (l in 1:N){
            alpha_old <- Alpha[l,]
            alpha_star <- Alpha_star[l,]
            y <- Y[l,]
            r_mh <- r_mh + lgamma(sum(alpha_star))- lgamma(sum(alpha_star) + sum(y))- lgamma(sum(alpha_old))+ lgamma(sum(alpha_old) + sum(y))
            r_mh <- r_mh + lgamma(alpha_star[j] + y[j]) - lgamma(alpha_star[j]) - lgamma(alpha_old[j] + y[j]) + lgamma(alpha_old[j])
          }
          # r_mh <- 0.5*log(2*pi) - a_beta*log(b_beta) + lgamma(a_beta) - lgamma(a_beta + 0.5) + (a_beta + 0.5)*log(b_beta + beta_old^2/2)
          sum_delta <- sum(Delta_temp[,j])
          # r_mh <- r_mh + lgamma(a_p + sum_delta) + lgamma(R-sum_delta + b_p) - lgamma(a_p + sum_delta + 1) - lgamma(R - sum_delta - 1 + b_p)
          # r_mh <- r_mh - lgamma(a_p+1) - lgamma(b_p) + lgamma(a_p) + lgamma(b_p+1)
          lb_star <- lbeta(a_p + sum_delta, R-sum_delta + b_p) - lbeta(a_p, b_p)
          lb_old <- lbeta(a_p + sum_delta + 1, R - sum_delta - 1 + b_p) - lbeta(a_p, b_p)
          r_mh <- r_mh + log(1-exp(lb_star))-lb_old
          r_mh <- r_mh - dnorm(beta_old, 0, 10, log = TRUE)
          r_mh <- r_mh + dnorm(beta_old,0, 1, log = TRUE)
          # print(r_mh)
          if (r_mh > log(runif(1))){
            Beta_temp[r,j] <- 0
            log_alpha <- log_alpha_star
            Alpha <- Alpha_star
          }else{
            Delta_temp[r,j] <- 1
          }
        }
      }
    }
    
    # Update Beta
    for (j in 1:P){
      for (r in 1:R){
        if (Delta_temp[r,j] == 1){
          beta_old <- Beta_temp[r,j]
          beta_star <- rnorm(1, mean = beta_old, sd = 0.5)
          # beta_star <- rtruncnorm(1,  mean = beta_old, sd = 0.5, low = -5, high = 5)
          # Beta_star <- Beta_temp
          # Beta_star[r,j] <- beta_star
          # Alpha_star <- exp(X_design %*% rbind(Beta_0_temp, Beta_star))
          # b <- Beta_temp[,j]
          # b[r] <- beta_star
          # a <- X_design %*% c(Beta_0_temp[j],b)
          # Alpha_star <- Alpha
          # Alpha_star[,j] <- exp(a)
          
          log_alpha_star <- log_alpha
          for (i in 1:N){
            log_alpha_star[i,j] <- log_alpha_star[i,j] - beta_old*X[i,r] + beta_star*X[i,r]
          }
          Alpha_star <- exp(log_alpha_star)
          
          r_mh <- 0
          for (l in 1:N){
            alpha_old <- Alpha[l,]
            alpha_star <- Alpha_star[l,]
            y <- Y[l,]
            r_mh <- r_mh + lgamma(sum(alpha_star))- lgamma(sum(alpha_star) + sum(y))- lgamma(sum(alpha_old))+ lgamma(sum(alpha_old) + sum(y))
            r_mh <- r_mh + lgamma(alpha_star[j] + y[j]) - lgamma(alpha_star[j]) - lgamma(alpha_old[j] + y[j]) + lgamma(alpha_old[j])
          }
          r_mh <- r_mh + (a_beta+0.5)*(log(b_beta + beta_old^2/2) - log(b_beta + beta_star^2/2))
          # r_mh <- r_mh + (a_beta+0.5)*(log(b_beta*2 + beta_old^2) - log(b_beta*2 + beta_star^2))
          if (r_mh > log(runif(1))){
            Beta_temp[r,j] <- beta_star
            log_alpha <- log_alpha_star
            Alpha <- Alpha_star
          }
        }
      }
    }
    
    # Update Beta_0
    for (j in 1:P){
      beta_0_old <- Beta_0_temp[j]
      beta_0_star <- rnorm(1, mean = beta_0_old, sd = tau)
      # Beta_0_star <- Beta_0_temp
      # Beta_0_star[j] <- beta_0_star
      # Alpha_star <- exp(X_design %*% rbind(Beta_0_star, Beta_temp))
      # a <- X_design %*% c(beta_0_star,Beta_temp[,j])
      # Alpha_star <- Alpha
      # Alpha_star[,j] <- exp(a)
      
      log_alpha_star <- log_alpha
      for (i in 1:N){
        log_alpha_star[i,j] <- log_alpha_star[i,j] - beta_0_old + beta_0_star
      }
      Alpha_star <- exp(log_alpha_star)
      
      r_mh <- 0
      for (l in 1:N){
        alpha_old <- Alpha[l,]
        alpha_star <- Alpha_star[l,]
        y <- Y[l,]
        r_mh <- r_mh + lgamma(sum(alpha_star))- lgamma(sum(alpha_star) + sum(y))- lgamma(sum(alpha_old))+ lgamma(sum(alpha_old) + sum(y))
        r_mh <- r_mh + lgamma(alpha_star[j] + y[j]) - lgamma(alpha_star[j]) - lgamma(alpha_old[j] + y[j]) + lgamma(alpha_old[j])
      }
      r_mh <- r_mh + (beta_0_old^2 - beta_0_star^2)/ 20
      if (r_mh > log(runif(1))){
        Beta_0_temp[j] <- beta_0_star
        log_alpha <- log_alpha_star
        Alpha <- Alpha_star
        acc <- acc + 1
      }
    }
    
    Beta_0_store[t,] <- Beta_0_temp
    Delta_store[t,,] <- Delta_temp
    Beta_store[t,,] <- Beta_temp
  }
  
  out <- NULL
  out$Beta_0_store <- Beta_0_store
  out$Delta_store <- Delta_store
  out$Beta_store <- Beta_store
  out$accept_rate <- acc/(T*P)
  
  return(out)
  
}


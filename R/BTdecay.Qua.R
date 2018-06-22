BTdecay.Qua <- function(dataframe, ability, theta, penalty.Qua, Lagrangian, decay.rate = 0, fixed = 1, iter = 100) {
  
  
  df <- dataframe
  u <- decay.rate
  v <- penalty.Qua
  n1 <- nrow(df)
  n <- nrow(ability) - 1
  
  fn <- function(ability) {
    
    s <- 0
    for(i in 1:(n-1)){
      for(j in (i+1):n){
        s1 <- theta[i, j] - ability[i] + ability[j]
        s <- s + Lagrangian[i, j] * s1 + v/2 * s1^2
      }
    }
    
    for (i in 1:n1) {
      a1 <- df[i, 1]
      a2 <- df[i, 2]
      C <- exp(-u * df[i, 5])
      x <- ability[a1]
      y <- ability[a2]
      p <- ability[n + 1] + x - y
      q <- exp(p)
      s <- s - (df[i, 3] * p - (df[i, 3] + df[i, 4]) * log(q + 1)) * C
    }
    s
  }
  
  gr <- function(ability) {
    
    Grad <- rep(0, n + 1)
    for(i in 1:n){
      Grad[i] <- v * (- sum(ability[-(n + 1)]) + n * ability[i] +
                          sum(theta[, i]) - sum(theta[i, ])) +
        sum(Lagrangian[, i]) - sum(Lagrangian[i, ])
    }
    for(i in 1:n1){
      a1 <- df[i, 1]
      a2 <- df[i, 2]
      C <- exp(-u * df[i, 5])
      x <- ability[a1]
      y <- ability[a2]
      p <- ability[n + 1] + x - y
      q <- exp(p)
      A <- -(df[i, 3] * (1/(q + 1)) + df[i, 4] * (-q/(q + 1))) * C
      Grad[a1] <- Grad[a1] + A
      Grad[a2] <- Grad[a2] - A
      Grad[n + 1] <- Grad[n + 1] + A
    }
    Grad
  }
  
  xa <- optimr::optimr(rep(0, n + 1), fn, gr = gr, method = "L-BFGS-B", control = list(maxit = iter))
  
  if(xa$convergence == 1){
    stop("Iterations diverge, please provide a smaller decay rate or more data")
  }
  
  ability[, 1] <- xa$par - xa$par[fixed]
  ability[n + 1, 1] <- xa$par[n + 1]
  ability
}


BTdecayLasso.step2 <- function(dataframe, ability, lambda, weight, decay.rate = 0, fixed = 1, thersh = 1e-5, iter = 100) {
  u <- decay.rate
  n <- nrow(ability) - 1
  theta <- matrix(0, nrow = n, ncol = n) 
  Lagrangian <- matrix(0, nrow = n, ncol = n) 
  ability[, 1] <- 0
  con <- matrix(NA, nrow = 0, ncol = 4)
  
  stop <- 0
  v <- 1
  j <- 1
  while (stop==0) {
    ability <- BTdecayLasso.step1(dataframe, ability, weight, Lagrangian, theta, v, lambda, 
                                  decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
    theta <- BTtheta(ability, weight, Lagrangian, v, lambda)
    Lagrangian0 <- BTLagrangian(Lagrangian, ability, theta, v)
    k <- sum(abs(Lagrangian0 - Lagrangian))
    if (k < thersh) {
      stop <- 1
    } else {
      Lagrangian <- Lagrangian0
      v <- max(Lagrangian^2)
    }
    s <- penaltyAmount(ability, weight)
    j <- j + 1
    con <- rbind(con, matrix(c(k, s, v, j), nrow = 1))
  }
  cat(s,'\n')
  
  ability0 <- ability
  ability0[, 1] <- 0
  BT <- BTdecay(dataframe, ability0, decay.rate = decay.rate, fixed = fixed, iter = iter)
  ability0 <- BT$ability
  s0 <- penaltyAmount(ability0, weight)
  p <- s/s0
  
  degree <- round(ability[1:n, 1], -log10(thersh)-1)
  degree <- length(unique(degree))
  
  output <- list(ability = ability, df = degree, penalty = p, convergence = con)
  output
}

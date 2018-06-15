#' @export

BTdecayLasso <- function(dataframe, ability, lambda = NULL, path = 89, decay.rate = 0, fixed = 1, thersh = 1e-5, max = 100, iter = 100) {
  u <- decay.rate
  n <- nrow(ability) - 1
  theta <- matrix(0, nrow = n, ncol = n) 
  Lagrangian <- matrix(0, nrow = n, ncol = n) 
  ability[, 1] <- 0
  Lambda <- seq(-0.5, -6, -0.0625)
  Lambda <- exp(Lambda)
  v <- 10
  weight <- BTLasso.weight(dataframe, ability, decay.rate = decay.rate, fixed = fixed, thersh = thersh, max = max, iter = iter)
  ability0 <- ability[, -1]
  l <- c()
  p <- c()
  
  BT <- BTdecay(dataframe, ability1, decay.rate = decay.rate, fixed = fixed, iter = iter)
  ability1 <- BT$ability
  s1 <- penaltyAmount(ability1, weight)
  l1 <- BTLikelihood(dataframe, ability1, decay.rate = decay.rate)
  df <- c()
  
  for (i in 1:length(Lambda)) {
    stop <- 0
    j <- 1
    
    while (stop==0) {
      ability <- BTdecayLasso.step1(dataframe, ability, weight, Lagrangian, theta, v, Lambda[i], 
                                    decay.rate = decay.rate, fixed = fixed, thersh = thersh, iter = iter)
      theta <- BTtheta(ability, weight, Lagrangian, v, Lambda[i])
      Lagrangian0 <- BTLagrangian(Lagrangian, ability, theta, v)
      k <- sum(abs(Lagrangian0 - Lagrangian))
      if (k < thersh) {
        stop <- 1
      } else {
        Lagrangian <- Lagrangian0
        v <- max(Lagrangian^2)
      }
      s0 <- penaltyAmount(ability, weight)
    }
    
    cat(s0, '\n')
    p0 <- s0/s1
    ability0 <- cbind(ability0, ability)
    l0 <- BTLikelihood(dataframe, BT$ability, decay.rate = decay.rate)
    l <- c(l, l0)
    p <- c(p, p0)
    
    degree <- round(ability[1:n, 1], -log10(thersh)-1)
    degree <- length(unique(degree))
    df <- c(df, degree)
  }
  
  ability0 <- cbind(ability0, ability1)
  l <- c(l, l1)
  p <- c(p, p1)
  df <- c(df, n)
  
  output <- list(ability = ability, likelihood = l, penalty = p, df = df, Lamda = c(Lambda, 0))
  output
}
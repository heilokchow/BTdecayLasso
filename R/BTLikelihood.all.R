BTLikelihood.all <- function(dataframe, ability, theta, penalty.Qua, weight, Lagrangian, lambda, decay.rate = 0){
  df <- dataframe
  v <- penalty.Qua
  u <- decay.rate
  n1 <- nrow(df)
  n <- nrow(ability) - 1
  s <- 0
  for(i in 1:n1){
    a1 <- df[i, 1]
    a2 <- df[i, 2]
    C <- exp(-u * df[i, 5])
    x <- ability[a1, 1]
    y <- ability[a2, 1]
    p <- ability[n + 1, 1] + x - y
    q <- exp(p)
    s <- s - (df[i, 3] * (p - log(q + 1)) + df[i, 4] * (-log(q + 1))) * C
  }
  theta1 <- theta
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      theta1[i, j] <- theta[i, j] - ability[i, 1] + ability[j, 1]
    }
  }
  s <- s + v/2 * sum(theta1^2) + sum(Lagrangian * theta1)
  s <- s + lambda * sum(abs(theta) * weight)
  s
}
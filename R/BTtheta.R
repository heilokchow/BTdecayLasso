BTtheta <- function(ability, weight, Lagrangian, penalty.Qua, lambda){
  n <- nrow(ability) - 1
  theta <- matrix(0, nrow = n, ncol = n) 
  v <- penalty.Qua
  
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      theta0 <- ability[i, 1] - ability[j, 1] - Lagrangian[i, j]/v
      theta[i, j] <- sign(theta0) * max(abs(theta0) - lambda * weight[i, j]/v, 0)
    }
  }
  return(theta)
}